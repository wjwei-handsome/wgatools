use crate::parser::{
    common::{AlignRecord, RecStat},
    maf::MAFReader,
    paf::PAFReader,
};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::{
    collections::HashMap,
    error::Error,
    io::{Read, Write},
};

/// Pair of reference and query as KEY
#[derive(Debug, Hash, PartialEq, Eq, Clone, Serialize, Deserialize, Default)]
struct Pair {
    ref_name: String,
    query_name: String,
    ref_size: u64,
    query_size: u64,
}

/// Statistic of a pair, Serialize for output
#[derive(Debug, Serialize, Deserialize, Default)]
struct Statistic {
    ref_name: String,
    query_name: String,
    ref_size: u64,
    query_size: u64,
    aligned_size: usize, // aggre by each block
    unaligned_size: u64,
    matched: usize,    // agg
    mismatched: usize, // agg
    ins_event: usize,  // agg
    del_event: usize,  // agg
    ins_size: usize,   // agg
    del_size: usize,   // agg
    identity: f32,
    similarity: f32,
    inv_event: usize,     // agg
    inv_size: f32,        // agg
    inv_ins_event: usize, // agg
    inv_ins_size: usize,  // agg
    inv_del_event: usize, // agg
    inv_del_size: usize,  // agg
}

// define a type for pair_stat
struct PairStat {
    pair: Pair,
    rec_stat: RecStat,
}

// stat for maf
pub fn stat_maf<R: Read + Send>(
    mut reader: MAFReader<R>,
    writer: &mut dyn Write,
    each: bool,
) -> Result<(), Box<dyn Error>> {
    let pair_stat_vec = reader
        .records()
        .par_bridge()
        .map(|rec| {
            let rec = match rec {
                Ok(r) => r,
                Err(e) => {
                    return Err(Box::new(e));
                }
            };
            Ok(stat_rec(&rec))
        })
        .flatten()
        .collect::<Vec<_>>();

    write_style_result(pair_stat_vec, writer, each)
}

// stat for paf
// TODO: impl the stat_rec for PAfRecord
pub fn stat_paf<R: Read + Send>(
    mut reader: PAFReader<R>,
    writer: &mut dyn Write,
    each: bool,
) -> Result<(), Box<dyn Error>> {
    let pair_stat_vec = reader
        .records()
        .par_bridge()
        .map(|rec| {
            let rec = match rec {
                Ok(r) => r,
                Err(e) => {
                    return Err(Box::new(e));
                }
            };
            Ok(stat_rec(&rec))
        })
        .flatten()
        .collect::<Vec<_>>();

    write_style_result(pair_stat_vec, writer, each)
}

fn write_style_result(
    pair_stat_vec: Vec<PairStat>,
    writer: &mut dyn Write,
    each: bool,
) -> Result<(), Box<dyn Error>> {
    let mut final_stat = match each {
        true => split_final(pair_stat_vec),
        false => merge_final_from_pair(pair_stat_vec),
    };
    final_stat.sort_by(|a, b| natord::compare(&a.ref_name, &b.ref_name));
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_writer(writer);
    for stat in final_stat {
        wtr.serialize(stat)?;
    }
    wtr.flush()?;
    Ok(())
}

// not merge
fn split_final(pair_stat_vec: Vec<PairStat>) -> Vec<Statistic> {
    // init final_stat
    let mut final_stat = Vec::new();
    for pair_stat in pair_stat_vec {
        let pair = pair_stat.pair;
        let rec_stat = pair_stat.rec_stat;
        // for rach record, init a Statistic
        let mut stat = Statistic {
            ref_name: pair.ref_name,
            query_name: pair.query_name,
            ref_size: pair.ref_size,
            query_size: pair.query_size,
            ..Default::default()
        };
        stat.aligned_size = rec_stat.aligned_size;
        stat.matched = rec_stat.matched;
        stat.mismatched = rec_stat.mismatched;
        stat.ins_event = rec_stat.ins_event;
        stat.del_event = rec_stat.del_event;
        stat.ins_size = rec_stat.ins_size;
        stat.del_size = rec_stat.del_size;
        stat.inv_ins_event = rec_stat.inv_ins_event;
        stat.inv_ins_size = rec_stat.inv_ins_size;
        stat.inv_del_event = rec_stat.inv_del_event;
        stat.inv_del_size = rec_stat.inv_del_size;
        stat.inv_event = rec_stat.inv_event;
        stat.inv_size = rec_stat.inv_size;
        stat.identity = stat.matched as f32 / stat.aligned_size as f32;
        stat.similarity = (stat.matched + stat.mismatched) as f32 / stat.aligned_size as f32;
        // push to final_stat
        final_stat.push(stat);
    }
    final_stat
}

// merge blocks in aggregation by ref_name
fn merge_final_from_pair(pair_stat_vec: Vec<PairStat>) -> Vec<Statistic> {
    // init final_stat
    let mut final_stat = Vec::new();

    let mut pair_stat_map = HashMap::new();
    for pair_stat in pair_stat_vec {
        let pair = pair_stat.pair;
        let rec_stat = pair_stat.rec_stat;
        pair_stat_map
            .entry(pair)
            .or_insert(Vec::new())
            .push(rec_stat);
    }
    for (pair, rec_stats) in pair_stat_map {
        // for rach record, init a Statistic
        let mut stat = Statistic {
            ref_name: pair.ref_name,
            query_name: pair.query_name,
            ref_size: pair.ref_size,
            query_size: pair.query_size,
            ..Default::default()
        };
        // aggregate by each record
        for rec_stat in rec_stats {
            stat.aligned_size += rec_stat.aligned_size;
            stat.matched += rec_stat.matched;
            stat.mismatched += rec_stat.mismatched;
            stat.ins_event += rec_stat.ins_event;
            stat.del_event += rec_stat.del_event;
            stat.ins_size += rec_stat.ins_size;
            stat.del_size += rec_stat.del_size;
            stat.inv_ins_event += rec_stat.inv_ins_event;
            stat.inv_ins_size += rec_stat.inv_ins_size;
            stat.inv_del_event += rec_stat.inv_del_event;
            stat.inv_del_size += rec_stat.inv_del_size;
            stat.inv_event += rec_stat.inv_event;
            stat.inv_size += rec_stat.inv_size;
        }
        // calculate the identity and similarity
        stat.unaligned_size = stat.ref_size - stat.aligned_size as u64;
        stat.identity = stat.matched as f32 / stat.aligned_size as f32;
        stat.similarity = (stat.matched + stat.mismatched) as f32 / stat.aligned_size as f32;
        // push to final_stat
        final_stat.push(stat);
    }
    final_stat
}

// stat a record to generate a PairStat
fn stat_rec<T: AlignRecord>(rec: &T) -> PairStat {
    // get pair
    let ref_name = rec.target_name();
    let query_name = rec.query_name();
    let ref_size = rec.target_length();
    let query_size = rec.query_length();
    let pair = Pair {
        ref_name: ref_name.to_string(),
        query_name: query_name.to_string(),
        ref_size,
        query_size,
    };

    // get rec_stat
    let rec_stat = rec.get_stat();

    PairStat { pair, rec_stat }
}
