use crate::{
    errors::WGAError,
    parser::{
        common::AlignRecord,
        paf::{PAFReader, PafRecord},
    },
};
use std::{
    collections::HashMap,
    io::{Read, Write},
};

// main function of generate pesudo MAF from PAF
pub fn trim_ovp<R: Read + Send>(
    mut reader: PAFReader<R>,
    writer: &mut dyn Write,
) -> Result<(), WGAError> {
    // 1. gourp by target
    let mut target_groupby_map: HashMap<String, Vec<PafRecord>> = HashMap::new();
    // TODO: could be parallel
    for rec in reader.records() {
        let rec = rec?;
        let target_name = rec.target_name().to_string();
        let rec_vec = target_groupby_map.entry(target_name).or_default();
        rec_vec.push(rec);
    }
    // TODO: could be parallel
    for (_, rec_vec) in target_groupby_map {
        trim_query(writer, rec_vec)?;
    }

    Ok(())
}

// filter for each query
fn trim_query(writer: &mut dyn Write, rec_vec: Vec<PafRecord>) -> Result<(), WGAError> {
    // groupby query name and sort by target start
    // [A,B,C,D1,D2,E] => {A:[A],B:[B],C:[C],D:[D1,D2],E:E}
    let mut query_groupby_map: HashMap<String, Vec<PafRecord>> = HashMap::new();
    // TODO: could be parallel
    for rec in rec_vec {
        let query_name = rec.query_name().to_string();
        let query_rec_vec = query_groupby_map.entry(query_name).or_default();
        // sort by query start
        let idx = query_rec_vec
            .binary_search_by(|probe| probe.target_start().cmp(&rec.target_start()))
            .unwrap_or_else(|e| e);
        query_rec_vec.insert(idx, rec);
    }
    let mut final_rec_vec: Vec<PafRecord> = Vec::new();
    // TODO: could be parallel
    for (_, rec_vec) in query_groupby_map {
        // remove overlap by target_start and target_end; keep the longest
        let mut rec_vec_iter = rec_vec.into_iter();
        let mut last_rec = rec_vec_iter.next().unwrap();
        for rec in rec_vec_iter {
            let align_size = rec.target_end() - rec.target_start();
            let last_align_size = last_rec.target_end() - last_rec.target_start();
            // remove overlap by target_start and target_end; keep the longest
            if rec.target_start() >= last_rec.target_end() {
                // no overlap
                final_rec_vec.push(last_rec);
                last_rec = rec;
            } else if align_size > last_align_size {
                // overlap, but longer
                last_rec = rec;
            } else {
                // overlap, but shorter
                continue;
            }
        }
        final_rec_vec.push(last_rec);
    }
    let mut pafwtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer);
    for rec in final_rec_vec {
        pafwtr.serialize(rec)?;
    }
    Ok(())
}
