// use wgalib::parser::common::FileFormat;
// use wgalib::parser::paf::PAFReader;
use wgalib::parser::maf::MAFReader;

fn main() {
    // let mut reader = PafReader::from_path("/Users/wjwei/Zm-CML333.paf").unwrap();
    //
    // reader.convert("test.chain", FileFormat::Chain);
    // let mut reader = PafReader::from_path("/Users/wjwei/Zm-CML333.paf").unwrap();
    // reader.convert("tes.blocks", FileFormat::Blocks);
    //
    let mut reader = MAFReader::from_path("data/output/test.maf").unwrap();
    let _header = &reader.header;
    for record in reader.records() {
        // let record = record.unwrap();
        // println!("{:?}", record);
        println!("{}", record.unwrap().get_cigar());
    }
    // println!("{:?}", reader.inner.buffer()); // buffer has been consumed
}
