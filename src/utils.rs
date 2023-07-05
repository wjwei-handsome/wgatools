use std::fs::File;
use std::io::{stdout, BufWriter, Write};

/// get a output writer including stdout and file writer
pub fn output_writer(outputpath: &str) -> Box<dyn Write> {
    if outputpath == "-" {
        Box::new(stdout())
    } else {
        let file = File::create(outputpath).unwrap();
        Box::new(BufWriter::new(file))
    }
}
