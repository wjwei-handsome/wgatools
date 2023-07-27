use crate::parser::chain::{ChainReader, ChainRecord};
use crate::parser::common::AlignRecord;
use crate::utils::output_writer;
use rust_htslib::faidx;
use std::io;

/// Convert a Chain Reader to output a MAF file
pub fn chain2maf<R: io::Read+Send>(
    chainreader: &mut ChainReader<R>,
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

    for chainrec in chainreader.records() {
        let chainrec = chainrec.unwrap();
        let t_name = chainrec.target_name();
        let t_start = chainrec.target_start() as usize;
        let t_end = (chainrec.target_end() - 1) as usize;
        let t_strand = chainrec.target_strand();
        let t_alilen = chainrec.target_end() - chainrec.target_start();
        let t_size = chainrec.target_length() as usize;
        let q_name = chainrec.query_name();
        let q_start = chainrec.query_start() as usize;
        let q_end = (chainrec.query_end() - 1) as usize;
        let q_strand = chainrec.query_strand();
        let q_size = chainrec.query_length() as usize;
        let q_alilen = chainrec.query_end() - chainrec.query_start();

        // write a a-line
        writer
            .write_all(format!("a\tscore={}\n", chainrec.header.chain_id).as_bytes())
            .unwrap();
        // fetch the sequence
        let mut whole_t_seq = t_reader.fetch_seq_string(t_name, t_start, t_end).unwrap();
        let mut whole_q_seq = q_reader.fetch_seq_string(q_name, q_start, q_end).unwrap();
        // nom the cigar string and insert the `-` to sequence
        let () = parse_chain_to_insert(&chainrec, &mut whole_t_seq, &mut whole_q_seq).unwrap();
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

/// Parse the Chain Data Lines to insert the `-` to sequence
fn parse_chain_to_insert(
    rec: &ChainRecord,
    t_seq: &mut String,
    q_seq: &mut String,
) -> io::Result<()> {
    let mut current_offset = 0;
    for dataline in &rec.lines {
        let match_len = dataline.size - current_offset;
        let ins_len = dataline.target_diff;
        let del_len = dataline.query_diff;
        current_offset += match_len;
        match ins_len {
            0 => {},
            _ => {
                let ins_str = "-".repeat(ins_len as usize);
                t_seq.insert_str(current_offset as usize, &ins_str);
                current_offset += ins_len;
            }
        }
        match del_len {
            0 => {},
            _ => {
                let del_str = "-".repeat(del_len as usize);
                q_seq.insert_str(current_offset as usize, &del_str);
                current_offset += del_len;
            }
        }
    }
    Ok(())
}