use crate::parser::cigar::cigar_cat_ext;
use crate::parser::common::{AlignRecord, Strand};
use crate::parser::maf::{MAFReader, MAFRecord};
use itertools::Itertools;
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
use rayon::prelude::*;
use std::io::{Read, Write};

use super::index::MafIndex;

///
///A example:
///
/// ACGATGCTAGCT---ACG
/// AC--TGCTAACTGGGACG
/// within alignment: snp | ins | del | tandem expansion | tandem contraction | Repeat expansion | Repeat contraction
/// between alignment: INS | DEL | Repeat expansion | Repeat contraction

fn build_header(sample_name: &str) -> Header {
    let svlen_id = infokey::SV_LENGTHS;
    let svlen_info = Map::<Info>::from(&svlen_id);

    let svtype_id = infokey::SV_TYPE;
    let svtype_info = Map::<Info>::from(&svtype_id);

    let end_id = infokey::END_POSITION;
    let end_info = Map::<Info>::from(&end_id);

    let inv_nest_id = "INV_NEST".parse::<infokey::Key>().unwrap();
    let inv_nest_info = Map::<Info>::new(
        Number::Count(1),
        infotype::String,
        "Varations nested within inversion",
    );

    let queryinfo_id = "QI".parse::<gtkey::Key>().unwrap();
    let queryinfo_info =
        Map::<Format>::new(Number::Count(1), fmttype::String, "Query informations");

    let gt_id = gtkey::GENOTYPE;
    let gt_format = Map::<Format>::from(&gt_id);

    Header::builder()
        .add_info(svlen_id, svlen_info)
        .add_info(svtype_id, svtype_info)
        .add_info(end_id, end_info)
        .add_info(inv_nest_id, inv_nest_info)
        .add_format(queryinfo_id, queryinfo_info)
        .add_format(gt_id, gt_format)
        .add_sample_name(sample_name)
        .build()
}

pub fn call_var_maf<R: Read + Send>(
    mafreader: &mut MAFReader<R>,
    mafindex: Option<MafIndex>,
    writer: &mut dyn Write,
    if_snp: bool,
    svlen_cutoff: u64,
    _between: bool,
    sample: Option<&str>,
) {
    let mut vcf_wtr = vcf::Writer::new(writer);
    let sample = sample.unwrap_or("sample");
    let mut header = build_header(sample);

    let mut mafrecords = mafreader
        .records()
        .map(|rec| rec.unwrap())
        .collect::<Vec<_>>();
    // if sort
    mafrecords.sort();
    let within_var_recs = mafrecords
        .par_iter()
        .map(|mafrec| call_within_var(mafrec, if_snp, svlen_cutoff))
        .flatten()
        .collect::<Vec<_>>();
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
            header
                .contigs_mut()
                .insert(name.parse().unwrap(), contigmap);
        }
    }

    vcf_wtr.write_header(&header).expect("TODO: panic message");
    for rec in within_var_recs {
        vcf_wtr
            .write_record(&header, &rec)
            .expect("TODO: panic message");
    }
}

fn get_variant_rec(
    chro: &str,
    pos: usize,
    ref_base: &str,
    alt_base: &str,
    info: Option<&str>,
    format: Option<&str>,
) -> Record {
    // let keys = "GT".parse().unwrap();
    // let genotypes = Genotypes::new(keys, vec![vec![Some(Value::from("1|1"))]]);

    let genotypes = match format {
        Some(format) => format.parse().unwrap(),
        None => "GT\t1|1".parse().unwrap(),
    };

    let infos: recinfo = match info {
        Some(info) => match info.parse() {
            Ok(infos) => infos,
            Err(_) => recinfo::default(),
        },
        None => recinfo::default(),
    };
    Record::builder()
        .set_chromosome(chro.parse().unwrap())
        .set_position(Position::from(pos))
        .set_reference_bases(ref_base.parse().unwrap())
        .set_alternate_bases(alt_base.parse().unwrap())
        .set_info(infos)
        .set_genotypes(genotypes)
        .build()
        .unwrap()
}

fn call_within_var(mafrec: &MAFRecord, if_snp: bool, svlen_cutoff: u64) -> Vec<Record> {
    // target:ACG-TTTGATGCTAGCT---ACG
    // query :ACCATTT--TGCTAACTGGGACG

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
        var_recs.push(record);
    }

    let t_seq_iter = mafrec.target_seq().chars();
    let q_seq_iter = mafrec.query_seq().chars();
    let group_by_iter = t_seq_iter
        .zip(q_seq_iter)
        .group_by(|(c1, c2)| cigar_cat_ext(c1, c2));

    let mut init_info = String::new();
    if strand == Strand::Negative {
        init_info.push_str("INV_NEST=TRUE;");
    }

    for (k, g) in group_by_iter.into_iter() {
        let len = g.count() as u64;
        match k {
            '=' => {
                target_current_offset += len;
                query_current_offset += len;
            }
            'I' => {
                if len > svlen_cutoff {
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
                    var_recs.push(record);
                }
                query_current_offset += len;
            }
            'D' => {
                if len > svlen_cutoff {
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
                    var_recs.push(record);
                }
                target_current_offset += len;
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
                        let record = get_variant_rec(
                            chro,
                            target_current_offset as usize + 1,
                            ref_base,
                            alt_base,
                            None,
                            None,
                        );
                        var_recs.push(record);
                        target_current_offset += 1;
                        query_current_offset += 1;
                    }
                } else {
                    query_current_offset += len;
                    target_current_offset += len;
                }
            }
            _ => {}
        }
    }
    var_recs
}
