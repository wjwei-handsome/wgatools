use crate::parser::chain::Header;
use crate::parser::cigar::parse_cigar_to_chain;
use crate::parser::paf::PAFReader;
use crate::utils::output_writer;
use std::io;
use std::io::Write;

/// Convert a PAF Reader to output a Chain file
pub fn paf2chains<R: io::Read>(pafreader: &mut PAFReader<R>, outputpath: &str) {
    // init writer
    let mut writer = output_writer(outputpath);

    // iterate over records and give a self-increasing chain-id
    for (id, record) in pafreader.records().enumerate() {
        let record = record.unwrap(); // TODO: handle parse error

        // transform record to Chain Header
        let mut header = Header::from(&record);

        // set chain id
        header.chain_id = id;

        // write header without newline
        writer.write_all(format!("{}", header).as_bytes()).unwrap(); // TODO: handle IO error

        // nom the cigar string and write to file
        // TODO: handle IResult error
        let (_, _) = parse_cigar_to_chain(&record, &mut writer).unwrap();

        // additional newline for standard chain format
        writer.write_all(b"\n\n").unwrap();
    }
    writer.flush().unwrap(); // TODO: handle IO error
}
