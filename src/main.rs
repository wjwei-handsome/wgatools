use wgalib::parser::common::FileFormat;
// use wgalib::parser::paf::PAFReader;
use wgalib::parser::maf::MAFReader;

fn main() {
    // let mut reader = PAFReader::from_path("/Users/wjwei/Zm-CML333.paf").unwrap();
    //
    // reader.convert("test.chain", FileFormat::Chain);
    // let mut reader = PAFReader::from_path("/Users/wjwei/Zm-CML333.paf").unwrap();
    // reader.convert("tes.blocks", FileFormat::Blocks);

    let mut reader = MAFReader::from_path("data/output/test.maf").unwrap();
    let _header = &reader.header;
    // for record in reader.records() {
    //     let record = record.unwrap();
    //     // println!("{:?}", record);
    //     println!("{}\nquery_start:{}\tquery_end:{}\ntarget_start:{}\ttarget_end:{}",
    //              record.get_cigar_string(),
    //              record.query_start(),record.query_end(),
    //     record.target_start(),record.target_end());
    // }
    reader.convert("test.paf", FileFormat::Paf);
    // println!("{:?}", reader.inner.buffer()); // buffer has been consumed
}
