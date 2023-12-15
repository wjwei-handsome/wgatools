use crate::parser::{
    common::AlignRecord,
    maf::{MAFReader, MAFWriter},
    paf::PAFReader,
};
use std::{
    error::Error,
    io::{Read, Write},
};

// filter paf
pub fn filter_paf<R: Read + Send>(
    mut reader: PAFReader<R>,
    writer: &mut dyn Write,
    min_block_size: u64,
    min_query_size: u64,
) -> Result<(), Box<dyn Error>> {
    let mut pafwtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer);
    for rec in reader.records() {
        let rec = match rec {
            Ok(r) => r,
            Err(e) => {
                return Err(Box::new(e));
            }
        };
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
) -> Result<(), Box<dyn Error>> {
    // init a MAFWriter
    let mut mafwtr = MAFWriter::new(writer);
    // write header
    mafwtr.write_header("#maf version=1.6 filter=..".to_owned());
    for rec in reader.records() {
        let rec = match rec {
            Ok(r) => r,
            Err(e) => {
                return Err(Box::new(e));
            }
        };
        let rec = filter_alignrec(&rec, min_block_size, min_query_size)?;
        // just write the record
        if let Some(rec) = rec {
            mafwtr.write_record(rec);
        }
    }
    Ok(())
}

// filter record, return Option
fn filter_alignrec<T: AlignRecord>(
    rec: &T,
    min_block_size: u64,
    min_query_size: u64,
) -> Result<Option<&T>, Box<dyn Error>> {
    let query_length = rec.query_length();
    let block_length = rec.target_align_size();

    if (block_length < min_block_size) | (query_length < min_query_size) {
        return Ok(None);
    }

    Ok(Some(rec))
}
