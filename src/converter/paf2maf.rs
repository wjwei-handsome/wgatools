use std::io;
use rust_htslib::faidx;
use crate::parser::common::AlignRecord;
use crate::parser::paf::{PAFReader, PafRecord};

pub fn paf2maf<R: io::Read>(pafreader: &mut PAFReader<R>, outputpath: &str) {
    let mut t_reader = faidx::Reader::from_path("data/t.fa").unwrap();
    let mut q_reader = faidx::Reader::from_path("data/q.fa").unwrap();
    let mut pafreader = PAFReader::from_path("tiny.paf").unwrap();
    let pafrec = pafreader.records().next().unwrap().unwrap();
    let t_name = &pafrec.target_name;
    let t_start = pafrec.target_start as usize;
    let t_end = (pafrec.target_end -1) as usize;
    let q_name = &pafrec.query_name;
    let q_start = pafrec.query_start as usize;
    let q_end = (pafrec.query_end -1) as usize;
    println!("t_region: {}:{}-{}", t_name, t_start, t_end);
    println!("q_region: {}:{}-{}", q_name, q_start, q_end);
    let whole_t_seq = t_reader.fetch_seq_string(t_name, t_start, t_end).unwrap();
    let whole_q_seq = q_reader.fetch_seq_string(q_name, q_start, q_end).unwrap();
    println!("t_seq:{}\nq_seq:{}", whole_t_seq,whole_q_seq);

    let cigar_string = pafrec.get_cigar_string()
}
