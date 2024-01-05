use rayon::iter::{ParallelBridge, ParallelIterator};

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
    // let mut cov_map: HashMap<String, Vec<usize>> = HashMap::new();
    // for rec in reader.records() {
    //     let rec = rec?;
    //     let target_name = rec.target_name().to_string();
    //     let target_length = rec.target_length() as usize;
    //     let cov_vec = cov_map.entry(target_name).or_insert(vec![0; target_length]);
    //     let cigar = rec.get_cigar_str()?;
    //     let start = rec.target_start() as usize;
    //     update_cov_vec(cov_vec, cigar, start)?;
    // }

    // parallel
    let cov_map = reader
        .records()
        .par_bridge()
        .try_fold(HashMap::new, |mut acc: HashMap<String, Vec<usize>>, rec| {
            let rec = rec?;
            let target_name = rec.target_name().to_string();
            let target_length = rec.target_length() as usize;
            let cov_vec = acc.entry(target_name).or_insert(vec![0; target_length]);
            let cigar = rec.get_cigar_str()?;
            let start = rec.target_start() as usize;
            update_cov_vec(cov_vec, cigar, start)?;
            Ok::<HashMap<String, Vec<usize>>, WGAError>(acc)
        })
        .try_reduce(HashMap::new, |mut acc, mut map| {
            // acc.extend(map.drain());
            // accumulate count
            // SHIT!!: WHY NOT JUST ITERATOR
            for (target, cov_vec) in map.drain() {
                let acc_vec = acc.entry(target).or_insert(vec![0; cov_vec.len()]);
                for (acc, cov) in acc_vec.iter_mut().zip(cov_vec) {
                    *acc += cov;
                }
            }
            Ok(acc)
        })?;

    // Output in BED format
    for (target, coverage) in cov_map {
        for (pos, count) in coverage.iter().enumerate() {
            writeln!(writer, "{}\t{}\t{}\t{}", target, pos, pos + 1, count)?
        }
    }

    // Output collapsed coverage
    // for (target, coverage) in cov_map {
    //     let mut last_count = 0;
    //     let mut last_pos = 0;
    //     for (pos, count) in coverage.iter().enumerate() {
    //         if *count != last_count {
    //             writeln!(writer, "{}\t{}\t{}\t{}", target, last_pos, pos, last_count)?;
    //             last_count = *count;
    //             last_pos = pos;
    //         }
    //     }
    //     writeln!(
    //         writer,
    //         "{}\t{}\t{}\t{}",
    //         target,
    //         last_pos,
    //         coverage.len(),
    //         last_count
    //     )?;
    // }
    Ok(())
}
