use crate::{
    errors::WGAError,
    parser::{cigar::update_cov_vec, common::AlignRecord, paf::PAFReader},
};
use std::{
    collections::HashMap,
    io::{Read, Write},
};

// main function of PAF Coverage
pub fn pafcov<R: Read + Send>(
    mut reader: PAFReader<R>,
    writer: &mut dyn Write,
) -> Result<(), WGAError> {
    let mut cov_map = HashMap::new();
    for rec in reader.records() {
        let rec = rec?;
        let target_name = rec.target_name().to_string();
        let target_length = rec.target_length() as usize;
        let cov_vec = cov_map.entry(target_name).or_insert(vec![0; target_length]);
        let cigar = rec.get_cigar_str()?;
        let start = rec.target_start() as usize;
        update_cov_vec(cov_vec, cigar, start)?;
    }
    // Output in BED format
    for (target, coverage) in cov_map {
        for (pos, count) in coverage.iter().enumerate() {
            writeln!(writer, "{}\t{}\t{}\t{}", target, pos, pos + 1, count)?
        }
    }
    Ok(())
}
