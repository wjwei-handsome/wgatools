use wgalib::parser::common::FileFormat;
use wgalib::parser::paf::PAFReader;
use wgalib::parser::maf::MAFReader;

fn main() {
    // let mut reader = PafReader::from_path("/Users/wjwei/Zm-CML333.paf").unwrap();
    //
    // reader.convert("test.chain", FileFormat::Chain);
    // let mut reader = PafReader::from_path("/Users/wjwei/Zm-CML333.paf").unwrap();
    // reader.convert("tes.blocks", FileFormat::Blocks);
    //
    let mut reader = MafReader::from_path("/Users/wjwei/Zm-CML333.maf").unwrap();
}
