use crate::parser::maf::MAFReader;
use crate::utils::output_writer;
use std::io;

/// Convert a MAF Reader to output a PAF file
pub fn maf2paf<R: io::Read>(mafreader: &mut MAFReader<R>, outputpath: &str) {
    // init writer and csv writer for deserializing
    let writer = output_writer(outputpath);
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);

    // iterate over records
    for record in mafreader.records() {
        let record = record.unwrap(); // TODO: handle parse error

        // nom the cigar string and write to file
    }
}