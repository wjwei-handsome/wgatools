use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::parser::{common::AlignRecord, maf::MAFReader};

#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub struct Pair {
    ref_name: String,
    query_name: String,
    ref_size: u64,
    query_size: u64,
}

#[derive(Debug, Serialize, Deserialize, Default)]
pub struct Statistic {
    ref_name: String,
    query_name: String,
    ref_size: u64,
    query_size: u64,
    aligned_size: u64, // aggre by each block
    unaligned_size: u64,
    matched: u64,    // agg
    mismatched: u64, // agg
    ins_event: u64,  // agg
    del_event: u64,  // agg
    ins_size: u64,   // agg
    del_size: u64,   // agg
    identity: f32,
    similarity: f32,
    inv_event: u64,
    inv_size: f32,
    inv_ins_event: u64, // agg
    inv_ins_size: u64,  // agg
    inv_del_event: u64, // agg
    inv_del_size: u64,  // agg
}

#[derive(Debug, Serialize, Deserialize, Default)]
pub struct RecStat {
    aligned_size: u64,  // aggre by each block
    matched: u64,       // agg
    mismatched: u64,    // agg
    ins_event: u64,     // agg
    del_event: u64,     // agg
    ins_size: u64,      // agg
    del_size: u64,      // agg
    inv_ins_event: u64, // agg
    inv_ins_size: u64,  // agg
    inv_del_event: u64, // agg
    inv_del_size: u64,  // agg
}

pub type PairStat = HashMap<Pair, Vec<RecStat>>;

fn gen_final_from_pair(pair_stat: &PairStat) -> Vec<Statistic> {
    let mut finals = Vec::new();
    for (pair, rec_stats) in pair_stat {
        let mut stat = Statistic::default();
        stat.ref_name = pair.ref_name.clone();
        stat.query_name = pair.query_name.clone();
        stat.ref_size = pair.ref_size;
        stat.query_size = pair.query_size;
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
        }
        stat.unaligned_size = stat.ref_size - stat.aligned_size;
        stat.identity = stat.matched as f32 / stat.aligned_size as f32;
        stat.similarity = (stat.matched + stat.mismatched) as f32 / stat.aligned_size as f32;
        stat.inv_size = stat.inv_ins_size as f32 / stat.inv_del_size as f32;
        finals.push(stat);
    }
    finals
}

pub fn stat() {
    let mut mafrdr = MAFReader::from_path("data/output/test_CML69.maf").unwrap();
    let mut pair_stat = PairStat::new();
    for rec in mafrdr.records() {
        let rec = rec.unwrap();
        stat_rec(&rec, &mut pair_stat);
    }
    let final_stat = gen_final_from_pair(&pair_stat);
    let stdout = std::io::stdout();
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_writer(stdout);
    final_stat.into_iter().for_each(|s| {
        wtr.serialize(s).unwrap(); // TODO: handle serialize error
    });
    wtr.flush().unwrap();
}

fn stat_rec<T: AlignRecord>(rec: &T, pair_stat: &mut PairStat) {
    // let mut pair_stat = PairStat::new();
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
    let mut stat = RecStat::default();
    stat.aligned_size = rec.target_align_size();
    stat.matched = 1;
    stat.mismatched = 1;
    stat.ins_event = 1;
    stat.del_event = 1;
    stat.ins_size = 1;
    stat.del_size = 1;
    stat.inv_ins_event = 1;
    stat.inv_ins_size = 1;
    stat.inv_del_event = 1;
    stat.inv_del_size = 1;

    // insert pair
    pair_stat.entry(pair.clone()).or_insert(Vec::new());
    pair_stat.get_mut(&pair).unwrap().push(stat);
}
