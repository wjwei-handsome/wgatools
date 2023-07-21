use wgalib::caller::test_call;
use wgalib::parser::common::FileFormat;
use wgalib::parser::maf::MAFReader;
// use wgalib::parser::paf::PAFReader;

fn main() {
    // let mut reader = MAFReader::from_path("data/test.maf").unwrap();
    // let _header = &reader.header;
    // reader.convert("tiny.bam", FileFormat::Bam);
    // let mut reader = MAFReader::from_path("/Users/wjwei/NAM2B73v5.Zm-CML333__TO__Zm-B73v5.maf").unwrap();
    // let _header = &reader.header;
    // reader.convert("test-out.paf", FileFormat::Paf);
    // let mut reader = MAFReader::from_path("data/test.maf").unwrap();
    // let _header = &reader.header;
    // reader.convert("tiny.chain", FileFormat::Chain);

    // let mut reader = PAFReader::from_path("/Users/wjwei/Zm-CML333.paf").unwrap();
    // reader.convert("test.chain", FileFormat::Chain);
    // let mut reader = PAFReader::from_path("/Users/wjwei/Zm-CML333.paf").unwrap();
    // reader.convert("tes.blocks", FileFormat::Blocks);
    // let mut reader = PAFReader::from_path("test.paf").unwrap();
    // reader.convert(
    //     "test-out.maf",
    //     FileFormat::Maf,
    //     Some("/Users/wjwei/Zm-B73v5.genome.fa"),
    //     Some("/Users/wjwei/Zm-CML333.genome.fa"),
    // );

    test_call();
}
