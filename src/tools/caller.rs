use crate::errors::WGAError;
use crate::parser::cigar::cst2cu;
use crate::parser::cigar::{cigar_cat_ext_caller, parse_cigar_str_tuple};
use crate::parser::common::{AlignRecord, Strand};
use crate::parser::maf::{MAFReader, MAFRecord, MAFSLine};
use crate::parser::paf::{PAFReader, PafRecord};
use crate::tools::index::MafIndex;
use itertools::Itertools;
use log::info;
use nom::bytes::complete::tag;
use nom::multi::fold_many1;
use noodles::vcf;
use noodles::vcf::{
    header::{
        record::value::{
            map::{format::Type as fmttype, info::Type as infotype, Contig, Format, Info},
            Map,
        },
        Number,
    },
    record::{
        genotypes::keys::key as gtkey, info::field::key as infokey, Info as recinfo, Position,
    },
    Header, Record,
};
use rust_htslib::faidx;
use rust_htslib::faidx::Reader as FaReader;
use std::io::{Read, Write};

// A example:
//
// ACGATGCTAGCT---ACG
// AC--TGCTAACTGGGACG
// within alignment: snp | ins | del | tandem expansion | tandem contraction | Repeat expansion | Repeat contraction
// between alignment: INS | DEL | Repeat expansion | Repeat contraction

// main function, it return a Result<(), WGAErr>
// NOTE: but other functions took anyhow, bucause noodles::vcf's error' organization is too complex
// and it will not be error in 99.9% cases
#[allow(clippy::too_many_arguments)]
pub fn call_var_maf<R: Read + Send>(
    mafreader: &mut MAFReader<R>,
    mafindex: Option<MafIndex>,
    writer: &mut dyn Write,
    if_snp: bool,
    svlen_cutoff: u64,
    _between: bool,
    sample: Option<&str>,
    query_name: Option<&str>,
    chunk_size: Option<usize>,
) -> Result<(), WGAError> {
    let mut vcf_wtr = vcf::Writer::new(writer);
    let sample = sample.unwrap_or("sample");
    let mut header = build_header(sample)?;
    // add contig to header
    add_header_contig(mafindex, &mut header)?;
    vcf_wtr.write_header(&header)?;

    for maf_result in mafreader.records() {
        let maf_record = maf_result?;
        let base_chunk_size = chunk_size.unwrap_or(1000000);
        info!(
            "Processing record: {} with {} chunk size ",
            maf_record.target_name(),
            base_chunk_size
        );
        let total_size = maf_record.slines[0].seq.len();

        let mut chunk_start = 0;
        let mut chunk_count = 0;
        while chunk_start < total_size {
            // find save chunk boundary for avoid large SVs
            let (safe_end, next_start) = find_safe_chunk_boundary(
                &maf_record,
                chunk_start,
                base_chunk_size,
                svlen_cutoff,
                total_size,
            )?;

            let actual_chunk_size = safe_end - chunk_start;
            chunk_count += 1;
            info!(
                "Processed chunk {}: start={}, end={}, size={}, progress={:.2}%",
                chunk_count,
                chunk_start,
                safe_end,
                actual_chunk_size,
                (safe_end as f64 / total_size as f64) * 100.0
            );

            // create a new MAFRecord for this chunk
            let mut chunk_record = create_chunk_record(&maf_record, chunk_start, safe_end)?;

            // call variants within this chunk and write to VCF
            let var_recs = call_within_var(&mut chunk_record, if_snp, svlen_cutoff, query_name)?;
            for rec in var_recs {
                vcf_wtr.write_record(&header, &rec)?;
            }

            chunk_start = next_start;
        }
        info!(
            "Finished processing record: {} chunks processed",
            chunk_count
        );
    }

    Ok(())
}

fn find_safe_chunk_boundary(
    record: &MAFRecord,
    start: usize,
    chunk_size: usize,
    svlen_cutoff: u64,
    total_size: usize,
) -> Result<(usize, usize), WGAError> {
    let proposed_end = (start + chunk_size).min(total_size);
    let mut current_gap_size = 0;
    let mut in_sv = false;
    let mut sv_start = 0;
    let mut safe_end = proposed_end;

    // check if the chunk contains a large structural variant
    for (pos, (ref_char, query_char)) in record.target_seq()[start..proposed_end]
        .chars()
        .zip(record.query_seq()[start..proposed_end].chars())
        .enumerate()
    {
        let abs_pos = start + pos;

        // check if this is a gap
        if ref_char == '-' || query_char == '-' {
            if !in_sv {
                in_sv = true;
                sv_start = abs_pos;
            }
            current_gap_size += 1;
        } else if in_sv {
            // if the gap is large enough, we consider it a structural variant
            if current_gap_size >= svlen_cutoff as usize {
                // if the gap is at the beginning of the chunk, we need to scan further to find the end of the SV
                if sv_start >= start {
                    // keep scanning until we find a non-gap character
                    safe_end = abs_pos;
                }
            }
            in_sv = false;
            current_gap_size = 0;
        }
    }

    // if the chunk ends with a structural variant, we need to scan further to find the end of the SV
    if in_sv && current_gap_size >= svlen_cutoff as usize {
        let mut end_pos = proposed_end;
        for (pos, (ref_char, query_char)) in record.target_seq()[proposed_end..]
            .chars()
            .zip(record.query_seq()[proposed_end..].chars())
            .enumerate()
        {
            if ref_char != '-' && query_char != '-' {
                end_pos = proposed_end + pos;
                break;
            }
        }
        safe_end = end_pos;
    }

    // return the safe end position and the next start position
    Ok((safe_end, safe_end))
}

fn create_chunk_record(
    original: &MAFRecord,
    start: usize,
    end: usize,
) -> Result<MAFRecord, WGAError> {
    let mut chunk = MAFRecord {
        score: original.score,
        slines: Vec::with_capacity(original.slines.len()),
        query_idx: original.query_idx,
    };

    for sline in &original.slines {
        let seq = &sline.seq[start..end];

        // get the new start position and alignment size
        let mut new_start = sline.start;
        let mut new_align_size = 0;

        // cal the bases
        for c in sline.seq[..start].chars() {
            if c != '-' {
                new_start += 1;
            }
        }

        // cal the alignment size
        for c in seq.chars() {
            if c != '-' {
                new_align_size += 1;
            }
        }

        chunk.slines.push(MAFSLine {
            mode: sline.mode,
            name: sline.name.clone(),
            start: new_start,
            align_size: new_align_size,
            strand: sline.strand,
            size: sline.size,
            seq: seq.to_string(),
        });
    }

    Ok(chunk)
}

#[allow(clippy::too_many_arguments)]
pub fn call_var_paf<R: Read + Send>(
    pafreader: &mut PAFReader<R>,
    t_fa_path: &str,
    q_fa_path: &str,
    writer: &mut dyn Write,
    if_snp: bool,
    svlen_cutoff: u64,
    _between: bool,
    sample: Option<&str>,
) -> Result<(), WGAError> {
    let mut vcf_wtr = vcf::Writer::new(writer);
    let sample = sample.unwrap_or("sample");
    let mut header = build_header(sample)?;

    // Process records sequentially to avoid sharing FASTA readers between threads
    let mut within_var_recs = Vec::new();
    let t_reader = faidx::Reader::from_path(t_fa_path)?;
    let q_reader = faidx::Reader::from_path(q_fa_path)?;

    // Collect and process records sequentially
    for pafrec_result in pafreader.records() {
        let pafrec = pafrec_result?;
        let var_recs = call_within_var_paf(&pafrec, if_snp, svlen_cutoff, &t_reader, &q_reader)?;
        within_var_recs.extend(var_recs);
    }

    // write VCF
    add_header_contig(None, &mut header)?;
    vcf_wtr.write_header(&header)?;
    for rec in within_var_recs {
        vcf_wtr.write_record(&header, &rec)?;
    }

    Ok(())
}

fn build_header(sample_name: &str) -> anyhow::Result<Header> {
    let svlen_id = infokey::SV_LENGTHS;
    let svlen_info = Map::<Info>::from(&svlen_id);

    let svtype_id = infokey::SV_TYPE;
    let svtype_info = Map::<Info>::from(&svtype_id);

    let end_id = infokey::END_POSITION;
    let end_info = Map::<Info>::from(&end_id);

    let inv_nest_id = "INV_NEST".parse::<infokey::Key>()?;
    let inv_nest_info = Map::<Info>::new(
        Number::Count(1),
        infotype::String,
        "Varations nested within inversion",
    );

    let queryinfo_id = "QI".parse::<gtkey::Key>()?;
    let queryinfo_info =
        Map::<Format>::new(Number::Count(1), fmttype::String, "Query informations");

    let gt_id = gtkey::GENOTYPE;
    let gt_format = Map::<Format>::from(&gt_id);

    Ok(Header::builder()
        .add_info(svlen_id, svlen_info)
        .add_info(svtype_id, svtype_info)
        .add_info(end_id, end_info)
        .add_info(inv_nest_id, inv_nest_info)
        .add_format(queryinfo_id, queryinfo_info)
        .add_format(gt_id, gt_format)
        .add_sample_name(sample_name)
        .build())
}

fn add_header_contig(mafindex: Option<MafIndex>, header: &mut Header) -> anyhow::Result<()> {
    if let Some(mafindex) = mafindex {
        let mut contig_vec: Vec<(String, u64)> = Vec::new();
        for (name, item) in mafindex {
            if item.ord == 0 {
                let size = item.size;
                contig_vec.push((name, size));
            }
            // natual sort by name use natord::compare
            contig_vec.sort_by(|a, b| natord::compare(&a.0, &b.0));
        }
        for (name, size) in contig_vec {
            let mut contigmap = Map::<Contig>::new();
            *contigmap.length_mut() = Some(size as usize);
            header.contigs_mut().insert(name.parse()?, contigmap);
        }
    }
    Ok(())
}

fn get_variant_rec(
    chro: &str,
    pos: usize,
    ref_base: &str,
    alt_base: &str,
    info: Option<&str>,
    format: Option<&str>,
) -> anyhow::Result<Record> {
    // let genotypes = Genotypes::new(keys, vec![vec![Some(Value::from("1|1"))]]);

    let genotypes = match format {
        Some(format) => format.parse()?,
        None => "GT\t1|1".parse()?,
    };

    let infos: recinfo = match info {
        Some(info) => info.parse().unwrap_or_default(),
        None => recinfo::default(),
    };
    Ok(Record::builder()
        .set_chromosome(chro.parse()?)
        .set_position(Position::from(pos))
        .set_reference_bases(ref_base.parse()?)
        .set_alternate_bases(alt_base.parse()?)
        .set_info(infos)
        .set_genotypes(genotypes)
        .build()?)
}

fn call_within_var(
    mafrec: &mut MAFRecord,
    if_snp: bool,
    svlen_cutoff: u64,
    query_name: Option<&str>,
) -> Result<Vec<Record>, WGAError> {
    // target:ACG-TTTGATGCTAGCT---ACG
    // query :ACCATTT--TGCTAACTGGGACG

    match query_name {
        Some(qname) => mafrec.set_query_idx_byname(qname)?,
        None => mafrec.set_query_idx(1),
    }

    let mut var_recs = Vec::new();

    let mut target_current_offset = mafrec.target_start();
    let mut query_current_offset = mafrec.query_start();

    let chro = mafrec.target_name();
    let q_chro = mafrec.query_name();
    let t_start = mafrec.target_start();
    let t_end = mafrec.target_end();
    let q_start = mafrec.query_start();
    let q_end = mafrec.query_end();

    let init_format: &str = "GT:QI\t1|1:";

    let mut t_seq_ref = mafrec.target_seq().to_string();
    t_seq_ref.retain(|c| c != '-');

    let mut q_seq_ref = mafrec.query_seq().to_string();
    q_seq_ref.retain(|c| c != '-');

    // add a inversion record if MafRecord's strand is '-'
    let strand = mafrec.query_strand();
    let format_surfix = match strand {
        Strand::Negative => 'N',
        Strand::Positive => 'P',
    };
    if strand == Strand::Negative {
        let ref_base = &t_seq_ref[0..1];
        let info = format!("SVTYPE=INV;END={}", t_end);
        let queryinfo = format!(
            "{}{}@{}@{}@{}",
            init_format, q_chro, q_start, q_end, format_surfix
        );
        let record = get_variant_rec(
            chro,
            target_current_offset as usize + 1,
            ref_base,
            "<INV>",
            // &id,
            Some(&info),
            Some(&queryinfo),
        );
        var_recs.push(record?);
    }

    let t_seq_iter = mafrec.target_seq().chars();
    let q_seq_iter = mafrec.query_seq().chars();
    let group_by_iter = t_seq_iter
        .zip(q_seq_iter)
        .group_by(|(c1, c2)| cigar_cat_ext_caller(c1, c2));

    let mut init_info = String::new();
    if strand == Strand::Negative {
        init_info.push_str("INV_NEST=TRUE;");
    }
    let mut after_m = false;
    for (k, g) in group_by_iter.into_iter() {
        let len = g.count() as u64;
        match k {
            '=' => {
                target_current_offset += len;
                query_current_offset += len;
                after_m = true;
            }
            'W' => {
                // do nothing
            }
            'I' => {
                if len > svlen_cutoff {
                    // This case for:
                    // t: ----A
                    // q: AAAAA
                    // This case for:
                    // t: TTGGG---
                    // q: TT---AAA
                    // This case for:
                    // t: GGG---
                    // q: ---AAA
                    if !after_m {
                        query_current_offset += len;
                        after_m = false;
                        continue;
                    }
                    let t_slice_start = (target_current_offset - t_start - 1) as usize;
                    let t_slice_end = t_slice_start + 1;

                    let q_slice_start = (query_current_offset - q_start - 1) as usize;
                    let q_slice_end = q_slice_start + len as usize + 1;

                    let info = format!(
                        "{}SVTYPE=INS;SVLEN={};END={}",
                        init_info, len, target_current_offset
                    );

                    let queryinfo = format!(
                        "{}{}@{}@{}@{}",
                        init_format,
                        q_chro,
                        query_current_offset,
                        query_current_offset + len,
                        format_surfix
                    );

                    let ref_base = &t_seq_ref[t_slice_start..t_slice_end];
                    let alt_base = &q_seq_ref[q_slice_start..q_slice_end];
                    let record = get_variant_rec(
                        chro,
                        target_current_offset as usize,
                        ref_base,
                        alt_base,
                        // &id,
                        Some(&info),
                        Some(&queryinfo),
                    );
                    var_recs.push(record?);
                }
                query_current_offset += len;
                after_m = false;
            }
            'D' => {
                if len > svlen_cutoff {
                    // for this case:
                    // t: AAAAA
                    // q: ----A
                    // for this case:
                    // t: TT---AAAAA
                    // q: TTGGG----A
                    // for this case:
                    // t: ---AAAAA
                    // q: GGG----A

                    // normal case:
                    // t: TTAAAAA
                    // q: TT----A
                    if !after_m {
                        target_current_offset += len;
                        after_m = false;
                        continue;
                    }

                    let t_slice_start = (target_current_offset - t_start - 1) as usize;
                    let t_slice_end = t_slice_start + len as usize + 1;

                    let q_slice_start = (query_current_offset - q_start - 1) as usize;
                    let q_slice_end = q_slice_start + 1;

                    let end = target_current_offset + len;
                    let info = format!("{}SVTYPE=DEL;SVLEN={};END={}", init_info, len, end);
                    let queryinfo = format!(
                        "{}{}@{}@{}@{}",
                        init_format,
                        q_chro,
                        query_current_offset,
                        query_current_offset,
                        format_surfix
                    );
                    // let id = format!("DEL{}", del_count);
                    let ref_base = &t_seq_ref[t_slice_start..t_slice_end];
                    let alt_base = &q_seq_ref[q_slice_start..q_slice_end];
                    let record = get_variant_rec(
                        chro,
                        target_current_offset as usize,
                        ref_base,
                        alt_base,
                        // &id,
                        Some(&info),
                        Some(&queryinfo),
                    );
                    var_recs.push(record?);
                }
                target_current_offset += len;
                after_m = false;
            }
            'X' => {
                if if_snp {
                    for _ in 0..len {
                        let t_slice_start = (target_current_offset - t_start) as usize;
                        let t_slice_end = t_slice_start + 1;

                        let q_slice_start = (query_current_offset - q_start) as usize;
                        let q_slice_end = q_slice_start + 1;

                        let ref_base = &t_seq_ref[t_slice_start..t_slice_end];
                        let alt_base = &q_seq_ref[q_slice_start..q_slice_end];

                        let queryinfo = format!(
                            "{}{}@{}@{}",
                            init_format, q_chro, query_current_offset, format_surfix
                        );
                        let record = get_variant_rec(
                            chro,
                            target_current_offset as usize + 1,
                            ref_base,
                            alt_base,
                            None,
                            Some(&queryinfo),
                        );
                        var_recs.push(record?);
                        target_current_offset += 1;
                        query_current_offset += 1;
                    }
                } else {
                    query_current_offset += len;
                    target_current_offset += len;
                }
                after_m = true;
            }
            _ => {}
        }
    }
    Ok(var_recs)
}

fn call_within_var_paf(
    pafrec: &PafRecord,
    if_snp: bool,
    svlen_cutoff: u64,
    t_fardr: &FaReader,
    q_fardr: &FaReader,
) -> Result<Vec<Record>, WGAError> {
    // init variant records
    let mut var_recs = Vec::new();

    let target_current_offset = pafrec.target_start();
    let query_current_offset = pafrec.query_start();

    let chro = pafrec.target_name();
    let q_chro = pafrec.query_name();
    let t_start = pafrec.target_start();
    let t_end = pafrec.target_end();
    let q_start = pafrec.query_start();
    let q_end = pafrec.query_end();

    let t_seq_string = pafrec.target_seq_with_fa(t_fardr)?;
    let q_seq_string = pafrec.query_seq_with_fa(q_fardr)?;

    let init_format: &str = "GT:QI\t1|1:";

    // add a inversion record if MafRecord's strand is '-'
    let strand = pafrec.query_strand();
    let format_surfix = match strand {
        Strand::Negative => 'N',
        Strand::Positive => 'P',
    };
    if strand == Strand::Negative {
        let ref_base = &t_seq_string[0..1];
        let info = format!("SVTYPE=INV;END={}", t_end);
        let queryinfo = format!(
            "{}{}@{}@{}@{}",
            init_format, q_chro, q_start, q_end, format_surfix
        );
        let record = get_variant_rec(
            chro,
            target_current_offset as usize + 1,
            ref_base,
            "<INV>",
            // &id,
            Some(&info),
            Some(&queryinfo),
        );
        var_recs.push(record?);
    }

    let cigar_string = pafrec.get_cigar_string()?;

    let (cigar, _tag) = tag("cg:Z:")(cigar_string.as_str())?;

    let mut t_pos = target_current_offset;
    let mut q_pos = query_current_offset;

    let mut init_info = String::new();
    if strand == Strand::Negative {
        init_info.push_str("INV_NEST=TRUE;");
    }
    let mut after_m = false;

    let (_, _) = fold_many1(parse_cigar_str_tuple, null, |res, cigarunit| {
        if res.is_ok() {
            // Parse CIGAR unit manually
            let cigarunit = cst2cu(cigarunit)?;
            let op = cigarunit.op;
            let len = cigarunit.len;
            match op {
                'M' | '=' => {
                    t_pos += len;
                    q_pos += len;
                    after_m = true;
                }
                'X' => {
                    if if_snp {
                        for _ in 0..len {
                            let t_slice_start = (t_pos - t_start) as usize;
                            let t_slice_end = t_slice_start + 1;

                            let q_slice_start = (q_pos - q_start) as usize;
                            let q_slice_end = q_slice_start + 1;

                            let ref_base = &t_seq_string[t_slice_start..t_slice_end];
                            let alt_base = &q_seq_string[q_slice_start..q_slice_end];

                            let queryinfo =
                                format!("{}{}@{}@{}", init_format, q_chro, q_pos, format_surfix);
                            let record = get_variant_rec(
                                chro,
                                t_pos as usize + 1,
                                ref_base,
                                alt_base,
                                None,
                                Some(&queryinfo),
                            );
                            var_recs.push(record?);
                            t_pos += 1;
                            q_pos += 1;
                        }
                    } else {
                        t_pos += len;
                        q_pos += len;
                    }
                    after_m = true;
                }
                'I' => {
                    if len > svlen_cutoff {
                        // This case for:
                        // t: ----A
                        // q: AAAAA
                        // This case for:
                        // t: TTGGG---
                        // q: TT---AAA
                        // This case for:
                        // t: GGG---
                        // q: ---AAA
                        if !after_m {
                            q_pos += len;
                            after_m = false;
                            return Ok(());
                        }
                        let t_slice_start = (t_pos - t_start - 1) as usize;
                        let t_slice_end = t_slice_start + 1;

                        let q_slice_start = (q_pos - q_start - 1) as usize;
                        let q_slice_end = q_slice_start + len as usize + 1;

                        let info = format!("{}SVTYPE=INS;SVLEN={};END={}", init_info, len, t_pos);

                        let queryinfo = format!(
                            "{}{}@{}@{}@{}",
                            init_format,
                            q_chro,
                            q_pos,
                            q_pos + len,
                            format_surfix
                        );

                        let ref_base = &t_seq_string[t_slice_start..t_slice_end];
                        let alt_base = &q_seq_string[q_slice_start..q_slice_end];
                        let record = get_variant_rec(
                            chro,
                            t_pos as usize,
                            ref_base,
                            alt_base,
                            Some(&info),
                            Some(&queryinfo),
                        );
                        var_recs.push(record?);
                    }
                    q_pos += len;
                    after_m = false;
                }
                'D' => {
                    if len > svlen_cutoff {
                        // for this case:
                        // t: AAAAA
                        // q: ----A
                        // for this case:
                        // t: TT---AAAAA
                        // q: TTGGG----A
                        // for this case:
                        // t: ---AAAAA
                        // q: GGG----A

                        // normal case:
                        // t: TTAAAAA
                        // q: TT----A
                        if !after_m {
                            t_pos += len;
                            after_m = false;
                            return Ok(());
                        }

                        let t_slice_start = (t_pos - t_start - 1) as usize;
                        let t_slice_end = t_slice_start + len as usize + 1;

                        let q_slice_start = (q_pos - q_start - 1) as usize;
                        let q_slice_end = q_slice_start + 1;

                        let info =
                            format!("{}SVTYPE=DEL;SVLEN={};END={}", init_info, len, t_pos + len);

                        let queryinfo = format!(
                            "{}{}@{}@{}@{}",
                            init_format, q_chro, q_pos, q_pos, format_surfix
                        );

                        let ref_base = &t_seq_string[t_slice_start..t_slice_end];
                        let alt_base = &q_seq_string[q_slice_start..q_slice_end];
                        let record = get_variant_rec(
                            chro,
                            t_pos as usize,
                            ref_base,
                            alt_base,
                            Some(&info),
                            Some(&queryinfo),
                        );
                        var_recs.push(record?);
                    }
                    t_pos += len;
                    after_m = false;
                }
                _ => return Err(WGAError::CigarOpInvalid(op.to_string())),
            }
        }
        res
    })(cigar)?;

    Ok(var_recs)
}

/// a phantom
fn null() -> Result<(), WGAError> {
    Ok(())
}
