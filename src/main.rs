use wgalib::parser::common::FileFormat;
use wgalib::parser::paf::PafReader;

fn main() {
    let mut reader = PafReader::from_path("/Users/wjwei/Zm-CML333.paf").unwrap();
    // for record in reader.records() {
    //     // let record = record.unwrap();
    //     println!("{:?}", record);
    // }
    reader.convert("test.chain", FileFormat::Chain);
    let mut reader = PafReader::from_path("/Users/wjwei/Zm-CML333.paf").unwrap();
    reader.convert("tes.blocks", FileFormat::Blocks);
    // match parse_cigar_to_alignment("cg:Z:20M5D15M".as_bytes()) {
    //     Ok((_, alignment)) => {
    //         println!("{:?}", alignment);
    //     }
    //     Err(e) => {
    //         println!("{:?}", e);
    //     }
    // }
    // paf2chain();
}
