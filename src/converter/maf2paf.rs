use crate::parser::common::AlignRecord;
use crate::parser::maf::MAFReader;
use crate::utils::output_writer;
use rayon::prelude::*;
use std::io;

/// Convert a MAF Reader to output a PAF file
pub fn maf2paf<R: io::Read + Send>(mafreader: &mut MAFReader<R>, outputpath: &str) {
    // init writer and csv writer for deserializing
    let writer = output_writer(outputpath);
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer);

    // multi-threading
    let pafrecords = mafreader
        .records()
        .par_bridge()
        .map(|record| {
            let mafrecord = record.unwrap(); // TODO: handle parse error
            mafrecord.convert2paf()
        })
        .collect::<Vec<_>>();
    // if we should sort pafrecords?
    pafrecords.into_iter().for_each(|record| {
        wtr.serialize(record).unwrap(); // TODO: handle serialize error
    });
    wtr.flush().unwrap(); // TODO: handle IO error
}
