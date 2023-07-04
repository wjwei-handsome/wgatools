use wgalib::parser::parse_cigar_to_alignment;

fn main() {
    // match PafReader::from_path("data/test.paf") {
    //     Ok(mut reader) => {
    //         for record in reader.records() {
    //             let record = record?;
    //             println!("{:?}", record);
    //         }
    //         // Ok(())
    //     }
    //     Err(e) => Err(e.into()),
    // };
    match parse_cigar_to_alignment("cg:Z:20M5D15M".as_bytes()) {
        Ok((_, alignment)) => {
            println!("{:?}", alignment);
        }
        Err(e) => {
            println!("{:?}", e);
        }
    }
}
