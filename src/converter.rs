// TODO: manage uses
use crate::errors::WGAError;
use crate::parser::chain::{ChainHeader, ChainReader, ChainRecord};
use crate::parser::cigar::{
    parse_cigar_to_blocks, parse_cigar_to_chain, parse_cigar_to_insert, parse_maf_seq_to_chain,
};
use crate::parser::common::{AlignRecord, Strand};
use crate::parser::maf::{MAFReader, MAFRecord, MAFSLine, MAFWriter};
use crate::parser::paf::PAFReader;
use crate::utils::reverse_complement;
use noodles::sam::header::record::value::map;
use noodles::sam::header::record::value::map::header::SortOrder;
use noodles::sam::record::ReadName;
use noodles::sam::{
    self as sam,
    header::record::value::{
        map::{Program, ReferenceSequence},
        Map,
    },
};
use rayon::prelude::*;
use rust_htslib::faidx;
use std::io::{Read, Write};
use std::num::NonZeroUsize;

/// Convert a MAF Reader to output a PAF file
pub fn maf2paf<R: Read + Send>(
    mafreader: &mut MAFReader<R>,
    writer: &mut dyn Write,
) -> Result<(), WGAError> {
    // init csv writer for deserializing
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer);

    // multi-threading
    let pafrecords = mafreader
        .records()
        .par_bridge()
        .map(|record| -> Result<_, WGAError> {
            let mafrecord = record?;
            mafrecord.convert2paf()
        })
        .collect::<Result<Vec<_>, WGAError>>()?;
    for pafrec in pafrecords {
        wtr.serialize(pafrec)?;
    }
    wtr.flush()?;
    Ok(())
}

/// Convert a MAF Reader to output a Chain file
pub fn maf2chain<R: Read + Send>(
    mafreader: &mut MAFReader<R>,
    writer: &mut Box<dyn Write>,
) -> Result<(), WGAError> {
    // iterate over records and give a self-increasing chain-id
    for (id, record) in mafreader.records().enumerate() {
        let record = record?;

        // transform record to Chain Header
        let mut header = ChainHeader::from(&record);

        // set chain id
        header.chain_id = id;

        // write header without newline
        writer.write_all(format!("{}", header).as_bytes())?;

        // nom the cigar string and write to file
        parse_maf_seq_to_chain(&record, writer)?;

        // additional newline for standard chain format
        writer.write_all(b"\n\n")?;
    }
    writer.flush()?;
    Ok(())
}

pub fn maf2sam<R: Read + Send>(
    _mafreader: &mut MAFReader<R>,
    writer: &mut Box<dyn Write>,
) -> Result<(), WGAError> {
    let mut sam_writer = sam::Writer::new(writer);
    let mut header = Map::<map::Header>::default();
    *header.sort_order_mut() = Some(SortOrder::Unsorted);
    let header = sam::Header::builder()
        .set_header(header)
        .add_reference_sequence(
            "sq0".parse()?,
            Map::<ReferenceSequence>::new(NonZeroUsize::try_from(8)?),
        )
        .add_reference_sequence(
            "sq1".parse()?,
            Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?),
        )
        .add_reference_sequence(
            "sq2".parse()?,
            Map::<ReferenceSequence>::new(NonZeroUsize::try_from(21)?),
        )
        .add_program("noodles-sam", Map::<Program>::default())
        .add_comment("an example SAM written by noodles-sam")
        .build();

    sam_writer.write_header(&header)?;
    let mut record = sam::alignment::Record::default();
    let read_name: ReadName = "sq2".parse()?;
    *record.read_name_mut() = Some(read_name);
    sam_writer.write_record(&header, &record)?;
    Ok(())
}

/// Convert a PAF Reader to output a Blocks file
pub fn paf2blocks<R: Read + Send>(
    pafreader: &mut PAFReader<R>,
    writer: &mut dyn Write,
) -> Result<(), WGAError> {
    // init writer and csv writer for deserializing
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);

    // iterate over records
    for record in pafreader.records() {
        let record = record?;
        // nom the cigar string and write to file
        parse_cigar_to_blocks(&record, &mut wtr)?;
    }
    wtr.flush()?;
    Ok(())
}

/// Convert a PAF Reader to output a Chain file
pub fn paf2chain<R: Read + Send>(
    pafreader: &mut PAFReader<R>,
    writer: &mut Box<dyn Write>,
) -> Result<(), WGAError> {
    // iterate over records and give a self-increasing chain-id
    for (id, record) in pafreader.records().enumerate() {
        let record = record?;

        // transform record to Chain Header
        let mut header = ChainHeader::from(&record);

        // set chain id
        header.chain_id = id;

        // write header without newline
        writer.write_all(format!("{}", header).as_bytes())?;

        // nom the cigar string and write to file
        parse_cigar_to_chain(&record, writer)?;

        // additional newline for standard chain format
        writer.write_all(b"\n\n")?;
    }
    writer.flush()?;
    Ok(())
}

/// Convert a PAF Reader to output a MAF file
pub fn paf2maf<R: Read + Send>(
    pafreader: &mut PAFReader<R>,
    writer: &mut dyn Write,
    t_fa_path: &str,
    q_fa_path: &str,
) -> Result<(), WGAError> {
    // get the target and query fasta reader
    let t_reader = faidx::Reader::from_path(t_fa_path)?;
    let q_reader = faidx::Reader::from_path(q_fa_path)?;

    // init a MAFWriter
    let mut mafwtr = MAFWriter::new(writer);

    // write header
    let header = format!(
        "#maf version=1.6 convert_from=paf t_seq_path={} q_seq_path={}",
        t_fa_path, q_fa_path
    );
    mafwtr.write_header(header)?;

    for pafrec in pafreader.records() {
        let pafrec = pafrec?;
        // get mapq as score
        let score = pafrec.mapq;
        // get target info
        let t_name = &pafrec.target_name;
        let t_start = pafrec.target_start;
        let t_end = pafrec.target_end - 1;
        let t_strand = pafrec.target_strand();
        let t_alilen = pafrec.target_end - pafrec.target_start;
        let t_size = pafrec.target_length;
        // get query info
        let q_name = &pafrec.query_name;
        let q_strand = pafrec.query_strand();
        let q_size = pafrec.query_length;
        let q_alilen = pafrec.query_end - pafrec.query_start;
        // NOTE: if negative strand, we should convert the start position
        let q_start = match q_strand {
            Strand::Positive => pafrec.query_start,
            Strand::Negative => q_size - pafrec.query_end,
        };

        // get seqs from indexed fasta files
        let mut whole_t_seq =
            t_reader.fetch_seq_string(t_name, t_start as usize, t_end as usize)?;
        let mut whole_q_seq = q_reader.fetch_seq_string(
            q_name,
            pafrec.query_start as usize,
            (pafrec.query_end - 1) as usize,
        )?;

        // reverse complement the query sequence if it is on the negative strand
        match q_strand {
            Strand::Positive => {}
            Strand::Negative => {
                whole_q_seq = reverse_complement(&whole_q_seq);
            }
        }
        // nom the cigar string and insert the `-` to sequence
        parse_cigar_to_insert(&pafrec, &mut whole_t_seq, &mut whole_q_seq)?;
        // get s-lines
        let t_sline = MAFSLine {
            mode: 's',
            name: t_name.to_string(),
            start: t_start,
            align_size: t_alilen,
            strand: t_strand,
            size: t_size,
            seq: whole_t_seq,
        };
        let q_sline = MAFSLine {
            mode: 's',
            name: q_name.to_string(),
            start: q_start,
            align_size: q_alilen,
            strand: q_strand,
            size: q_size,
            seq: whole_q_seq,
        };
        // get maf record
        let mafrec = MAFRecord {
            score,
            slines: vec![t_sline, q_sline],
        };
        // write maf record
        mafwtr.write_record(&mafrec)?;
    }
    Ok(())
}

/// Convert a Chain Reader to output a MAF file
pub fn chain2maf<R: Read + Send>(
    chainreader: &mut ChainReader<R>,
    writer: &mut dyn Write,
    t_fa_path: &str,
    q_fa_path: &str,
) -> Result<(), WGAError> {
    // get the target and query fasta reader
    let t_reader = faidx::Reader::from_path(t_fa_path)?;
    let q_reader = faidx::Reader::from_path(q_fa_path)?;

    // init a MAFWriter
    let mut mafwtr = MAFWriter::new(writer);

    // write header
    let header = format!(
        "#maf version=1.6 convert_from=chain t_seq_path={} q_seq_path={}",
        t_fa_path, q_fa_path
    );
    mafwtr.write_header(header)?;

    for chainrec in chainreader.records()? {
        let chainrec = chainrec?;
        // 255 as score
        let score = 255;
        // get target info
        let t_name = chainrec.target_name();
        let t_start = chainrec.target_start();
        let t_end = chainrec.target_end() - 1;
        let t_strand = chainrec.target_strand();
        let t_alilen = chainrec.target_end() - chainrec.target_start();
        let t_size = chainrec.target_length();
        // get query info
        let q_name = chainrec.query_name();
        let q_strand = chainrec.query_strand();
        let q_size = chainrec.query_length();
        let q_alilen = chainrec.query_end() - chainrec.query_start();
        // NOTE: if negative strand, we should convert the start position
        let q_start = match q_strand {
            Strand::Positive => chainrec.query_start(),
            Strand::Negative => q_size - chainrec.query_end(),
        };

        // get seqs from indexed fasta files
        let mut whole_t_seq =
            t_reader.fetch_seq_string(t_name, t_start as usize, t_end as usize)?;
        let mut whole_q_seq = q_reader.fetch_seq_string(
            q_name,
            chainrec.query_start() as usize,
            (chainrec.query_end() - 1) as usize,
        )?;

        // reverse complement the query sequence if it is on the negative strand
        match q_strand {
            Strand::Positive => {}
            Strand::Negative => {
                whole_q_seq = reverse_complement(&whole_q_seq);
            }
        }
        // read chain dataline and insert the `-` to sequence
        parse_chain_to_insert(&chainrec, &mut whole_t_seq, &mut whole_q_seq)?;
        // get s-lines
        let t_sline = MAFSLine {
            mode: 's',
            name: t_name.to_string(),
            start: t_start,
            align_size: t_alilen,
            strand: t_strand,
            size: t_size,
            seq: whole_t_seq,
        };
        let q_sline = MAFSLine {
            mode: 's',
            name: q_name.to_string(),
            start: q_start,
            align_size: q_alilen,
            strand: q_strand,
            size: q_size,
            seq: whole_q_seq,
        };
        // get maf record
        let mafrec = MAFRecord {
            score,
            slines: vec![t_sline, q_sline],
        };
        // write maf record
        mafwtr.write_record(&mafrec)?;
    }
    Ok(())
}

/// Parse the Chain Data Lines to insert the `-` to sequence
fn parse_chain_to_insert(
    rec: &ChainRecord,
    t_seq: &mut String,
    q_seq: &mut String,
) -> Result<(), WGAError> {
    let mut current_offset = 0;
    for dataline in &rec.lines {
        let ins_len = dataline.target_diff;
        let del_len = dataline.query_diff;
        current_offset += dataline.size;
        match ins_len {
            0 => {}
            _ => {
                let ins_str = "-".repeat(ins_len as usize);
                t_seq.insert_str(current_offset as usize, &ins_str);
                current_offset += ins_len;
            }
        }
        match del_len {
            0 => {}
            _ => {
                let del_str = "-".repeat(del_len as usize);
                q_seq.insert_str(current_offset as usize, &del_str);
                current_offset += del_len;
            }
        }
    }
    Ok(())
}

/// Convert a Chain Reader to output a PAF file
pub fn chain2paf<R: Read + Send>(
    chainreader: &mut ChainReader<R>,
    writer: &mut dyn Write,
) -> Result<(), WGAError> {
    // init csv writer for deserializing
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer);

    // multi-threading
    let pafrecords = chainreader
        .records()?
        .par_bridge()
        .map(|record| -> Result<_, WGAError> {
            let chainrecord = record?;
            chainrecord.convert2paf()
        })
        .collect::<Result<Vec<_>, WGAError>>()?;
    // if we should sort pafrecords?
    for pafrec in pafrecords {
        wtr.serialize(pafrec)?;
    }
    wtr.flush()?;
    Ok(())
}
