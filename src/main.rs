use wgalib::parser::common::FileFormat;
use wgalib::parser::paf::PAFReader;
// use wgalib::parser::maf::MAFReader;

fn main() {
    // let mut reader = PAFReader::from_path("/Users/wjwei/Zm-CML333.paf").unwrap();
    // reader.convert("test.chain", FileFormat::Chain);
    // let mut reader = PAFReader::from_path("/Users/wjwei/Zm-CML333.paf").unwrap();
    // reader.convert("tes.blocks", FileFormat::Blocks);
    let mut reader = PAFReader::from_path("tiny.paf").unwrap();
    reader.convert("tint-out.maf", FileFormat::Maf);

    // let mut reader = MAFReader::from_path("data/test.maf").unwrap();
    // let _header = &reader.header;
    // reader.convert("tiny.bam", FileFormat::Bam);
    // let mut reader = MAFReader::from_path("data/test.maf").unwrap();
    // let _header = &reader.header;
    // reader.convert("tiny.paf", FileFormat::Paf);
    // let mut reader = MAFReader::from_path("data/test.maf").unwrap();
    // let _header = &reader.header;
    // reader.convert("tiny.chain", FileFormat::Chain);
    // test_call();

}
