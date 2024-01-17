use crate::{
    errors::WGAError,
    parser::{
        cigar::gen_pesudo_maf_by_cigar,
        common::AlignRecord,
        paf::{PAFReader, PafRecord},
    },
};
use rust_htslib::faidx;
use std::{
    collections::HashMap,
    io::{Read, Write},
};

// main function of generate pesudo MAF from PAF
pub fn generate_pesudo_maf<R: Read + Send>(
    mut reader: PAFReader<R>,
    writer: &mut dyn Write,
    fa_path: &Option<String>,
    target: &Option<String>,
) -> Result<(), WGAError> {
    // 1. gourp by target and sort by target_start
    let mut final_map: HashMap<String, Vec<PafRecord>> = HashMap::new();
    // if Some(target), only output the target
    if let Some(target) = target {
        for rec in reader.records() {
            let rec = rec?;
            if rec.target_name() == target {
                let rec_vec = final_map.entry(target.to_string()).or_default();
                // sort by target start
                let idx = rec_vec
                    .binary_search_by(|probe| probe.target_start().cmp(&rec.target_start()))
                    .unwrap_or_else(|e| e);
                rec_vec.insert(idx, rec);
            }
        }
    } else {
        for rec in reader.records() {
            let rec = rec?;
            let target_name = rec.target_name().to_string();
            let rec_vec = final_map.entry(target_name).or_default();
            // sort by target start
            let idx = rec_vec
                .binary_search_by(|probe| probe.target_start().cmp(&rec.target_start()))
                .unwrap_or_else(|e| e);
            rec_vec.insert(idx, rec);
        }
    }
    // write header
    writeln!(writer, "##maf version=1")?;
    // 2. for each target, generate a pesudo maf
    for (target_name, rec_vec) in final_map {
        write_pmaf(writer, rec_vec, &target_name, fa_path)?;
    }

    Ok(())
}

/// RecRegions is a struct to store a Merged PafRecord and its query regions
struct RecRegions {
    rec: PafRecord,
    q_regions: Option<Vec<(u64, u64)>>,
}

// merge sorted rec_vec to a merged rec and query regions
fn merge_query_recs_map(
    query_recs_map: HashMap<String, Vec<PafRecord>>,
) -> Result<Vec<RecRegions>, WGAError> {
    let mut pafrec_vec = Vec::new();
    for (_, rec_vec) in query_recs_map {
        if rec_vec.len() == 1 {
            // push single rec without any change
            pafrec_vec.push(RecRegions {
                rec: rec_vec.into_iter().next().unwrap(),
                q_regions: None,
            });
            continue;
        } else {
            // merge the rec_vec to a single rec
            // new target_start is the min of target_start
            // new target_end is the max of target_end
            // new cigar is the merge of cigar, padding by `{n}D`, n is the distance between two rec

            // take first element as the new rec
            let mut q_regions = Vec::new();
            let mut rec_vec_iter = rec_vec.into_iter();
            let mut new_rec = match rec_vec_iter.next() {
                Some(rec) => rec,
                None => {
                    return Err(WGAError::Other(anyhow::anyhow!(
                        "rec_vec is empty, please contact the author"
                    )))
                }
            };
            let mut new_cigar = String::new();
            new_cigar.push_str(new_rec.get_cigar_str()?);
            let mut last_target_end = new_rec.target_end();
            // push first rec's query region
            let q_region = (new_rec.query_start(), new_rec.query_end());
            q_regions.push(q_region);
            for rec in rec_vec_iter {
                // push query region
                let q_region = (rec.query_start(), rec.query_end());
                q_regions.push(q_region);
                // trim "cg:Z:"
                let pure_cigar = rec.get_cigar_str()?.trim_start_matches("cg:Z:");
                // padding by `{n}D`
                let target_start = rec.target_start();
                let gap_len = target_start - last_target_end;
                new_cigar.push_str(&format!("{}D", gap_len));
                new_cigar.push_str(pure_cigar);
                // update last_target_end
                last_target_end = rec.target_end();
            }
            // update rec's cigar and target_end
            new_rec.tags = vec![new_cigar];
            new_rec.target_end = last_target_end;
            // push new merged rec
            pafrec_vec.push(RecRegions {
                rec: new_rec,
                q_regions: Some(q_regions),
            });
        }
    }
    Ok(pafrec_vec)
}

// write pseudo maf by merged rec and query regions
fn write_pmaf(
    writer: &mut dyn Write,
    rec_vec: Vec<PafRecord>,
    target_name: &str,
    fa_path: &Option<String>,
) -> Result<(), WGAError> {
    // groupby query name and sort by query start
    // [A,B,C,D1,D2,E] => {A:[A],B:[B],C:[C],D:[D1,D2],E:E}
    let mut query_recs_map: HashMap<String, Vec<PafRecord>> = HashMap::new();
    for rec in rec_vec {
        let query_name = rec.query_name().to_string();
        let query_rec_vec = query_recs_map.entry(query_name).or_default();
        // sort by query start
        let idx = query_rec_vec
            .binary_search_by(|probe| probe.query_start().cmp(&rec.query_start()))
            .unwrap_or_else(|e| e);
        query_rec_vec.insert(idx, rec);
    }

    // merge the query recs
    let merged_recregions_vec = merge_query_recs_map(query_recs_map)?;

    // start output
    writeln!(writer, "a score=0")?;

    // init first_flag to true
    let mut first_flag = true;
    // if fa_path is specified, use fa_path to fetch sequence
    // else use `N``,`0`,`1` to fill sequence
    let true_base = fa_path.is_some();

    for rec_regions in merged_recregions_vec {
        let rec = rec_regions.rec;

        // write target s-line
        let target_size = rec.target_length();
        let target_start = rec.target_start();
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

        // start write query s-line
        let query_name = rec.query_name();
        let query_size = rec.query_length();
        let cigar = rec.get_cigar_str()?;
        write!(
            writer,
            "s\t{}\t0\t{}\t+\t{}\t",
            query_name, query_size, query_size
        )?;
        // fill the head with '-'
        for _ in 0..target_start {
            write!(writer, "-")?;
        }

        // if merged rec, join query sequence
        let mut q_seq = match rec_regions.q_regions {
            // single rec
            None => {
                // get raw query sequence
                let q_start = rec.query_start();
                let q_end = rec.query_end();
                get_sline_seq(fa_path, query_name, (q_start, q_end), false)?
            }
            // merged rec
            Some(q_regions) => {
                // join query sequence
                let mut q_seq = String::new();
                for q_region in q_regions {
                    let q_seq_part = get_sline_seq(fa_path, query_name, q_region, false)?;
                    q_seq.push_str(&q_seq_part);
                }
                q_seq
            }
        };
        // modify query sequence by cigar
        gen_pesudo_maf_by_cigar(cigar, &mut q_seq, true_base)?;
        // write query sequence
        write!(writer, "{}", q_seq)?;
        // fill the tail with '-' if no secendary paf record
        let seq_line_len = q_seq.len() as u64;
        for _ in 0..(target_size - seq_line_len - target_start) {
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
