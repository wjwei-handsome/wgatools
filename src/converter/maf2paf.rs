use crate::parser::common::AlignRecord;
use crate::parser::maf::MAFReader;
use crate::utils::output_writer;
use rayon::prelude::*;
use std::io;

/// Convert a MAF Reader to output a PAF file
pub fn maf2paf<R: io::Read+Send>(mafreader: &mut MAFReader<R>, outputpath: &str) {
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
    let start0 = std::time::Instant::now();
    let pafrecords = mafreader
        .records()
        .par_bridge()
        .into_par_iter()
        .map(|record| {
            let mafrecord = record.unwrap(); // TODO: handle parse error
            mafrecord.convert2paf()
        })
        .collect::<Vec<_>>();
    println!("convert maf to paf: {:?}", start0.elapsed());
    // let start1 = std::time::Instant::now();
    // let mafrecords = mafreader
    //     .records()
    //     .map(|rec| rec.unwrap())
    //     .collect::<Vec<_>>();
    // println!("read maf records: {:?}", start1.elapsed());
    // let start2 = std::time::Instant::now();
    // let pafrecords = mafrecords
    //     .par_iter()
    //     .map(|record| record.convert2paf())
    //     .collect::<Vec<_>>();
    // println!("convert maf to paf: {:?}", start2.elapsed());
    let start3 = std::time::Instant::now();
    pafrecords.into_iter().for_each(|record| {
        wtr.serialize(record).unwrap(); // TODO: handle serialize error
    });
    wtr.flush().unwrap(); // TODO: handle IO error
    println!("write paf records: {:?}", start3.elapsed());
}
