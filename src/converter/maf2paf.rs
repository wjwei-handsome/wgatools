use crate::parser::common::AlignRecord;
use crate::parser::maf::MAFReader;
use crate::utils::output_writer;
use rayon::prelude::*;
use std::io;

/// Convert a MAF Reader to output a PAF file
pub fn maf2paf<R: io::Read>(mafreader: &mut MAFReader<R>, outputpath: &str) {
    // init writer and csv writer for deserializing
    let writer = output_writer(outputpath);
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer);

    // iterate over records
    // for record in mafreader.records() {
    //     let mafrecord = record.unwrap(); // TODO: handle parse error
    //
    //     let pafrecord = mafrecord.convert2paf();
    //     wtr.serialize(pafrecord).unwrap(); // TODO: handle serialize error
    // }
    // wtr.flush().unwrap(); // TODO: handle IO error

    // multi-threading
    let mafrecords = mafreader
        .records()
        .map(|rec| rec.unwrap())
        .collect::<Vec<_>>();
    let pafrecords = mafrecords
        .par_iter()
        .map(|record| record.convert2paf())
        .collect::<Vec<_>>();
    pafrecords.into_iter().for_each(|record| {
        wtr.serialize(record).unwrap(); // TODO: handle serialize error
    });
    wtr.flush().unwrap(); // TODO: handle IO error
}
