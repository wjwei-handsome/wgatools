use crate::{
    errors::WGAError,
    parser::{
        cigar::gen_pesudo_maf_by_cigar,
        common::{AlignRecord, Strand},
        paf::{PAFReader, PafRecord},
    },
    utils::reverse_complement,
};
use rayon::prelude::*;
use rust_htslib::faidx;
use std::{
    collections::HashMap,
    io::{Read, Write},
};

// main function of generate pesudo MAF from PAF
pub fn generate_pesudo_maf<R: Read + Send>(
    mut reader: PAFReader<R>,
    out_dir: &str,
    fa_path: &Option<String>,
    target: &Option<String>,
) -> Result<(), WGAError> {
    // 1. gourp by target
    let mut target_groupby_map: HashMap<String, Vec<PafRecord>> = HashMap::new();
    // if Some(target), only output the target
    if let Some(target) = target {
        for rec in reader.records() {
            let rec = rec?;
            if rec.target_name() == target {
                let rec_vec = target_groupby_map.entry(target.to_string()).or_default();
                rec_vec.push(rec);
            }
        }
    } else {
        for rec in reader.records() {
            let rec = rec?;
            let target_name = rec.target_name().to_string();
            let rec_vec = target_groupby_map.entry(target_name).or_default();
            rec_vec.push(rec);
        }
    }

    // output in parallel, each target will be output to a file
    // let mut handles = Vec::new();
    // for (target_name, rec_vec) in final_map {
    //     let mut out_path = out_dir.to_string();
    //     out_path.push_str(&format!("/{}.maf", target_name));
    //     let fa_path = fa_path.clone();
    //     let handle = std::thread::spawn(move || {
    //         let mut writer = std::fs::File::create(out_path)?;
    //         write_pmaf(&mut writer, rec_vec, &target_name, &fa_path)?;
    //         Ok::<(), WGAError>(())
    //     });
    //     handles.push(handle);
    // }
    // // wait for all thread to finish
    // for handle in handles {
    //     handle.join().unwrap()?;
    // }

    target_groupby_map
        .into_par_iter()
        .map(|(target_name, rec_vec)| {
            let mut out_path = out_dir.to_string();
            out_path.push_str(&format!("/{}.maf", target_name));
            let fa_path = fa_path.clone();
            let mut writer = std::fs::File::create(out_path)?;
            write_pmaf(&mut writer, rec_vec, &target_name, &fa_path)?;
            Ok::<(), WGAError>(())
        })
        .collect::<Result<Vec<_>, WGAError>>()?;

    Ok(())
}

// write pseudo maf by merged rec and query regions
fn write_pmaf(
    writer: &mut dyn Write,
    rec_vec: Vec<PafRecord>,
    target_name: &str,
    fa_path: &Option<String>,
) -> Result<(), WGAError> {
    // groupby query name and sort by target start
    // [A,B,C,D1,D2,E] => {A:[A],B:[B],C:[C],D:[D1,D2],E:E}
    let mut query_groupby_map: HashMap<String, Vec<PafRecord>> = HashMap::new();
    for rec in rec_vec {
        let query_name = rec.query_name().to_string();
        let query_rec_vec = query_groupby_map.entry(query_name).or_default();
        // sort by query start
        let idx = query_rec_vec
            .binary_search_by(|probe| probe.target_start().cmp(&rec.target_start()))
            .unwrap_or_else(|e| e);
        query_rec_vec.insert(idx, rec);
    }

    // start output
    writeln!(writer, "a score=0")?;

    // init first_flag to true
    let mut first_flag = true;
    // if fa_path is specified, use fa_path to fetch sequence
    // else use `N``,`0`,`1` to fill sequence
    let true_base = fa_path.is_some();
    // init first_query_rec_flag to true
    let mut target_size = 0;
    for (query_name, rec_vec) in query_groupby_map {
        let mut first_query_flag = true;
        let mut last_target_end = 0;
        for rec in rec_vec {
            // write target s-line
            target_size = rec.target_length();
            // let target_start = rec.target_start();
            if first_flag {
                write!(
                    writer,
                    "s\t{}\t0\t{}\t+\t{}\t",
                    target_name, target_size, target_size
                )?;
                let whole_t_seq = get_sline_seq(fa_path, target_name, (0, target_size), true)?;
                writeln!(writer, "{}", whole_t_seq)?;
                first_flag = false;
            }
            // start write query s-line if first_query_flag is true
            if first_query_flag {
                // start write query s-line
                let query_name = rec.query_name();
                let query_size = rec.query_length();
                write!(
                    writer,
                    "s\t{}\t0\t{}\t+\t{}\t",
                    query_name, query_size, query_size
                )?;
                // fill the head with '-'
                // for _ in 0..target_start {
                //     write!(writer, "-")?;
                // }
                // last_target_end = rec.target_end();
            }
            // process rest query recs
            let q_start = rec.query_start();
            let q_end = rec.query_end();

            // fill the gap between two query recs from second query rec
            // if !first_query_flag {
            let gap_len = rec.target_start() - last_target_end;
            for _ in 0..gap_len {
                write!(writer, "-")?;
            }
            last_target_end = rec.target_end();
            // }

            let mut q_seq = get_sline_seq(fa_path, &query_name, (q_start, q_end), false)?;
            // reverse complement the query sequence if it is on the negative strand
            match rec.query_strand() {
                Strand::Positive => {}
                Strand::Negative => {
                    q_seq = reverse_complement(&q_seq)?;
                }
            }
            // modify query sequence by cigar
            let cigar = rec.get_cigar_str()?;
            gen_pesudo_maf_by_cigar(cigar, &mut q_seq, true_base)?;
            // write modified query sequence
            write!(writer, "{}", q_seq)?;
            first_query_flag = false;
        }
        // fill the tail with '-'
        let tail_len = target_size - last_target_end;
        for _ in 0..tail_len {
            write!(writer, "-")?;
        }
        // new line
        writeln!(writer)?;
    }

    // final new line
    writeln!(writer)?;
    Ok(())
}

// get query sequence by query region if fa path is specified
// else return empty string
fn get_sline_seq(
    fa_path: &Option<String>,
    name: &str,
    region: (u64, u64),
    target: bool,
) -> Result<String, WGAError> {
    match fa_path {
        Some(path) => {
            let fa_reader = faidx::Reader::from_path(path)?;
            let q_start = region.0 as usize;
            let q_end = region.1 as usize - 1; // SHIT! FAIDX
            let raw_q_seq = fa_reader.fetch_seq_string(name, q_start, q_end)?;
            Ok(raw_q_seq)
        }
        None => {
            if target {
                // if target sequence, return `NNN..`
                Ok("N".repeat((region.1 - region.0) as usize))
            } else {
                Ok(String::new())
            }
        }
    }
}
