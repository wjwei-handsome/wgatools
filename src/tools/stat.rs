use crate::{
    errors::WGAError,
    parser::{
        common::{AlignRecord, RecStat},
        maf::MAFReader,
        paf::PAFReader,
    },
};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::{
    collections::HashMap,
    io::{Read, Write},
};

/// Pair of reference and query as KEY
#[derive(Debug, Hash, PartialEq, Eq, Clone, Serialize, Deserialize, Default)]
struct Pair {
    ref_name: String,
    ref_size: u64,
    query_name: String,
    query_size: u64,
}

/// Statistic of a pair, Serialize for output
#[derive(Debug, Serialize, Deserialize, Default)]
pub struct Statistic {
    pub ref_name: String,
    pub ref_size: u64,
    pub ref_start: u64,
    pub query_name: String,
    pub query_size: u64,
    pub query_start: u64,
    pub aligned_size: usize, // aggre by each block
    pub unaligned_size: u64,
    pub identity: f32,
    pub similarity: f32,
    pub matched: usize,       // agg
    pub mismatched: usize,    // agg
    pub ins_event: usize,     // agg
    pub del_event: usize,     // agg
    pub ins_size: usize,      // agg
    pub del_size: usize,      // agg
    pub inv_event: usize,     // agg
    pub inv_size: f32,        // agg
    pub inv_ins_event: usize, // agg
    pub inv_ins_size: usize,  // agg
    pub inv_del_event: usize, // agg
    pub inv_del_size: usize,  // agg
}

// define a type for pair_stat
struct PairStat {
    pair: Pair,
    ref_start: u64,
    query_start: u64,
    rec_stat: RecStat,
}

// stat for maf
pub fn stat_maf<R: Read + Send>(
    mut reader: MAFReader<R>,
    writer: &mut dyn Write,
    each: bool,
    query_name: Option<&str>,
) -> Result<(), WGAError> {
    let pair_stat_vec = reader
        .records()
        .par_bridge()
        .try_fold(Vec::new, |mut acc, result_rec| {
            let mut rec = result_rec?;
            if let Some(qname) = query_name {
                rec.set_query_idx_by_name(qname)?;
            }
            acc.push(stat_rec(&rec)?);
            Ok::<Vec<PairStat>, WGAError>(acc)
        })
        .try_reduce(Vec::new, |mut acc, mut vec| {
            acc.append(&mut vec);
            Ok(acc)
        })?;

    write_style_result(pair_stat_vec, writer, each)
}

// stat for paf
pub fn stat_paf<R: Read + Send>(
    mut reader: PAFReader<R>,
    writer: &mut dyn Write,
    each: bool,
) -> Result<(), WGAError> {
    let pair_stat_vec = reader
        .records()
        .par_bridge()
        .try_fold(Vec::new, |mut acc, rec| {
            acc.push(stat_rec(&rec?)?);
            Ok::<Vec<PairStat>, WGAError>(acc)
        })
        .try_reduce(Vec::new, |mut acc, mut vec| {
            acc.append(&mut vec);
            Ok(acc)
        })?;

    write_style_result(pair_stat_vec, writer, each)
}

fn write_style_result(
    pair_stat_vec: Vec<PairStat>,
    writer: &mut dyn Write,
    each: bool,
) -> Result<(), WGAError> {
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
            ref_size: pair.ref_size,
            ref_start: pair_stat.ref_start,
            query_name: pair.query_name,
            query_size: pair.query_size,
            query_start: pair_stat.query_start,
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
        pair_stat_map.entry(pair).or_insert(Vec::new()).push((
            rec_stat,
            pair_stat.ref_start,
            pair_stat.query_start,
        ));
    }
    for (pair, rec_stats) in pair_stat_map {
        // for rach record, init a Statistic
        let mut stat = Statistic {
            ref_name: pair.ref_name,
            ref_size: pair.ref_size,
            ref_start: pair.ref_size,
            query_name: pair.query_name,
            query_size: pair.query_size,
            query_start: pair.query_size,
            ..Default::default()
        };
        // aggregate by each record
        for (rec_stat, r_s, q_s) in rec_stats {
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
            // judge if ref_start and query_start is the smallest
            if r_s < stat.ref_start {
                stat.ref_start = r_s;
            }
            if q_s < stat.query_start {
                stat.query_start = q_s;
            }
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
fn stat_rec<T: AlignRecord>(rec: &T) -> Result<PairStat, WGAError> {
    // get pair
    let ref_name = rec.target_name();
    let ref_size = rec.target_length();
    let ref_start = rec.target_start();
    let query_name = rec.query_name();
    let query_size = rec.query_length();
    let query_start = rec.query_start();
    let pair = Pair {
        ref_name: ref_name.to_string(),
        ref_size,
        query_name: query_name.to_string(),
        query_size,
    };

    // get rec_stat
    let rec_stat = rec.get_stat()?;

    Ok(PairStat {
        pair,
        rec_stat,
        ref_start,
        query_start,
    })
}
