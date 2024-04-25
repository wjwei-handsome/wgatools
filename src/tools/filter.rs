use crate::{
    errors::WGAError,
    parser::{
        chain::ChainReader,
        common::AlignRecord,
        maf::{MAFReader, MAFWriter},
        paf::PAFReader,
    },
};
use rayon::prelude::*;
use std::{
    collections::HashMap,
    io::{Read, Write},
};

// filter chain
pub fn filter_chain<R: Read + Send>(
    mut reader: ChainReader<R>,
    writer: &mut dyn Write,
    min_block_size: u64,
    min_query_size: u64,
) -> Result<(), WGAError> {
    for rec in reader.records()? {
        let rec = rec?;
        let rec = filter_alignrec(&rec, min_block_size, min_query_size)?;
        // just write the record
        if let Some(rec) = rec {
            let chainheader = &rec.header;
            let chainblocks = &rec.lines;
            writer.write_all(format!("{}", chainheader).as_bytes())?;
            for dataline in chainblocks {
                writer.write_all(format!("{}", dataline).as_bytes())?;
            }
            // additional newline for standard chain format
            writer.write_all(b"\n\n")?;
        }
    }
    Ok(())
}

// filter paf
pub fn filter_paf<R: Read + Send>(
    mut reader: PAFReader<R>,
    writer: &mut dyn Write,
    min_block_size: u64,
    min_query_size: u64,
) -> Result<(), WGAError> {
    let mut pafwtr = csv::WriterBuilder::new()
        .flexible(true)
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer);
    for rec in reader.records() {
        let rec = rec?;
        let rec = filter_alignrec(&rec, min_block_size, min_query_size)?;
        // just write the record
        if let Some(rec) = rec {
            pafwtr.serialize(rec)?;
        }
    }
    Ok(())
}

// filter maf
pub fn filter_maf<R: Read + Send>(
    mut reader: MAFReader<R>,
    writer: &mut dyn Write,
    min_block_size: u64,
    min_query_size: u64,
) -> Result<(), WGAError> {
    // init a MAFWriter
    let mut mafwtr = MAFWriter::new(writer);
    // write header
    let header = format!(
        "#maf version=1.6 filter=blocksize>={} querysize>={}",
        min_block_size, min_query_size
    );
    mafwtr.write_header(header)?;
    for rec in reader.records() {
        let rec = rec?;
        let rec = filter_alignrec(&rec, min_block_size, min_query_size)?;
        // just write the record
        if let Some(rec) = rec {
            mafwtr.write_record(rec)?;
        }
    }
    Ok(())
}

// filter record, return Option
fn filter_alignrec<T: AlignRecord>(
    rec: &T,
    min_block_size: u64,
    min_query_size: u64,
) -> Result<Option<&T>, WGAError> {
    let query_length = rec.query_length();
    let block_length = rec.target_align_size();

    // if in condition, return None
    if (block_length < min_block_size) | (query_length < min_query_size) {
        return Ok(None);
    }

    Ok(Some(rec))
}

// main function of filter query-target pairs
pub fn filter_paf_align_pair<R: Read + Send>(
    mut reader: PAFReader<R>,
    writer: &mut dyn Write,
    filt_align_size: u64,
) -> Result<(), WGAError> {
    // parallel read and groupby
    let (align_size_sum_map, all_recs) = reader
        .records()
        .par_bridge()
        .try_fold(
            || (HashMap::new(), Vec::new()),
            |(mut align_size_sum_map, mut all_recs), rec| {
                let rec = rec?;
                let q_name = rec.query_name().to_string();
                let t_name = rec.target_name().to_string();
                let apx_align_size = rec.target_align_size();

                let key = (q_name, t_name);
                let entry = align_size_sum_map.entry(key).or_insert(0);
                *entry += apx_align_size;
                all_recs.push(rec);
                Ok::<_, WGAError>((align_size_sum_map, all_recs))
            },
        )
        .try_reduce(
            || (HashMap::new(), Vec::new()),
            |(mut align_size_sum_map1, mut all_recs1), (align_size_sum_map2, mut all_recs2)| {
                for (key, value) in align_size_sum_map2 {
                    let entry = align_size_sum_map1.entry(key).or_insert(0);
                    *entry += value;
                }
                all_recs1.append(&mut all_recs2);
                Ok((align_size_sum_map1, all_recs1))
            },
        )?;

    let mut pafwtr = csv::WriterBuilder::new()
        .flexible(true)
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer);
    // filter by align_size_sum
    for rec in all_recs {
        let q_name = rec.query_name().to_string();
        let t_name = rec.target_name().to_string();
        let key = (q_name, t_name);
        let align_size_sum = align_size_sum_map.get(&key).unwrap();
        if *align_size_sum >= filt_align_size {
            pafwtr.serialize(rec)?;
        }
    }
    Ok(())
}
