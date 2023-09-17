use crate::parser::common::AlignRecord;
use crate::parser::maf::MAFReader;


// fn region_intersect((s1: u64, e1: u64),(s2: u64, e2: u64)) {
//
// }

pub fn maf_extractor() {
    let input_chro = "Zm-B73v5.chr8";
    let input_start = 181121825;
    let input_end = 181131825;
    let input_length = input_start - input_end;
    let path = "data/output/test.maf";
    let mut reader = MAFReader::from_path(path).unwrap();

    for mafrec in reader.records() {
        let mafrec = mafrec.unwrap();
        let t_name = mafrec.target_name();
        let q_name = mafrec.query_name();
        if t_name == input_chro {
            let t_start = mafrec.target_start();
            let t_end = mafrec.target_end();
            if t_start <= input_start && t_end >= input_end {
                let cut_start_index = input_start - t_start;
                let cut_end_index = input_end - t_start;
                let t_seq = mafrec.target_seq();
                let q_seq = mafrec.query_seq();
                let t_seq = &t_seq[cut_start_index as usize..cut_end_index as usize];
                let q_seq = &q_seq[cut_start_index as usize..cut_end_index as usize];
                let new_q_start = mafrec.query_start() + cut_start_index;
                println!("s\t{}\t{}\t{}\t{}\t{}\t{}", t_name, input_start, input_length, mafrec.target_strand(), mafrec.target_length(), t_seq);
                println!("s\t{}\t{}\t{}\t{}\t{}\t{}", q_name, new_q_start, input_length, mafrec.query_strand(), mafrec.query_length(), q_seq);
            } else {
                todo!()
            }

        }
    }

}
