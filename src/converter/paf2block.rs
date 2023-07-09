use crate::parser::cigar::parse_cigar_to_blocks;
use crate::parser::paf::PAFReader;
use crate::utils::output_writer;
use std::io;

/// Convert a PAF Reader to output a Blocks file
pub fn paf2blocks<R: io::Read>(pafreader: &mut PAFReader<R>, outputpath: &str) {
    // init writer and csv writer for deserializing
    let writer = output_writer(outputpath);
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);

    // iterate over records
    for record in pafreader.records() {
        let record = record.unwrap(); // TODO: handle parse error

        // nom the cigar string and write to file
        // TODO: handle IResult error
        let (_, _) = parse_cigar_to_blocks(&record, &mut wtr).unwrap();
    }
    wtr.flush().unwrap(); // TODO: handle IO error
}
