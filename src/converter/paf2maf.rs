use crate::parser::cigar::parse_cigar_to_insert;
use crate::parser::common::AlignRecord;
use crate::parser::paf::PAFReader;
use crate::utils::output_writer;
use rust_htslib::faidx;
use std::io;
use std::io::Write;

pub fn paf2maf<R: io::Read>(
    pafreader: &mut PAFReader<R>,
    outputpath: &str,
    t_fa_path: Option<&str>,
    q_fa_path: Option<&str>,
) {
    // get the target and query fasta reader
    let t_fa_path = match t_fa_path {
        Some(path) => path,
        None => panic!("Please provide a target fasta file path"),
    };
    let q_fa_path = match q_fa_path {
        Some(path) => path,
        None => panic!("Please provide a query fasta file path"),
    };
    let t_reader = faidx::Reader::from_path(t_fa_path).unwrap();
    let q_reader = faidx::Reader::from_path(q_fa_path).unwrap();

    // init writer
    let mut writer = output_writer(outputpath);

    // write maf header
    writer.write_all(b"##maf\tversion=1\n").unwrap();

    for pafrec in pafreader.records() {
        let pafrec = pafrec.unwrap();
        let t_name = &pafrec.target_name;
        let t_start = pafrec.target_start as usize;
        let t_end = (pafrec.target_end - 1) as usize;
        let t_strand = pafrec.target_strand();
        let t_alilen = pafrec.target_end - pafrec.target_start;
        let t_size = pafrec.target_length as usize;
        let q_name = &pafrec.query_name;
        let q_start = pafrec.query_start as usize;
        let q_end = (pafrec.query_end - 1) as usize;
        let q_strand = pafrec.query_strand();
        let q_size = pafrec.query_length as usize;
        let q_alilen = pafrec.query_end - pafrec.query_start;

        // write a a-line
        writer
            .write_all(format!("a\tscore={}\n", pafrec.mapq).as_bytes())
            .unwrap();
        // fetch the sequence
        let mut whole_t_seq = t_reader.fetch_seq_string(t_name, t_start, t_end).unwrap();
        let mut whole_q_seq = q_reader.fetch_seq_string(q_name, q_start, q_end).unwrap();
        // nom the cigar string and insert the `-` to sequence
        let (_, _) = parse_cigar_to_insert(&pafrec, &mut whole_t_seq, &mut whole_q_seq).unwrap();
        // write a s-line
        writer
            .write_all(
                format!(
                    "s\t{}\t{}\t{}\t{}\t{}\t",
                    t_name, t_start, t_alilen, t_strand, t_size
                )
                .as_bytes(),
            )
            .unwrap();
        writer.write_all(whole_t_seq.as_bytes()).unwrap();
        writer.write_all(b"\n").unwrap();
        writer
            .write_all(
                format!(
                    "s\t{}\t{}\t{}\t{}\t{}\t",
                    q_name, q_start, q_alilen, q_strand, q_size
                )
                .as_bytes(),
            )
            .unwrap();
        writer.write_all(whole_q_seq.as_bytes()).unwrap();
        writer.write_all(b"\n").unwrap();
        writer.write_all(b"\n").unwrap();
    }
}
