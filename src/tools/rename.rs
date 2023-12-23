use crate::{
    errors::WGAError,
    parser::maf::{MAFReader, MAFWriter},
};
use std::io::{Read, Write};
// filter maf
pub fn rename_maf<R: Read + Send>(
    mut reader: MAFReader<R>,
    writer: &mut dyn Write,
    prefixs: Vec<&str>,
) -> Result<(), WGAError> {
    // init a MAFWriter
    let mut mafwtr = MAFWriter::new(writer);
    // write header
    let header = format!("#maf version=1.6 rename={}", prefixs.join(";"));
    mafwtr.write_header(header)?;
    for rec in reader.records() {
        let mut rec = rec?;
        rec.rename(&prefixs)?;
        mafwtr.write_record(&rec)?;
    }
    Ok(())
}
