use crate::parser::chain::ChainReader;
use crate::parser::common::AlignRecord;
use crate::utils::output_writer;
use rayon::prelude::*;
use std::io;

/// Convert a Chain Reader to output a PAF file
pub fn chain2paf<R: io::Read + Send>(chainreader: &mut ChainReader<R>, outputpath: &str) {
    // init writer and csv writer for deserializing
    let writer = output_writer(outputpath);
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer);

    // multi-threading
    let pafrecords = chainreader
        .records()
        .par_bridge()
        .map(|record| {
            let chainrecord = record.unwrap(); // TODO: handle parse error
            chainrecord.convert2paf()
        })
        .collect::<Vec<_>>();
    // if we should sort pafrecords?
    pafrecords.into_iter().for_each(|record| {
        wtr.serialize(record).unwrap(); // TODO: handle serialize error
    });
    wtr.flush().unwrap(); // TODO: handle IO error
}
