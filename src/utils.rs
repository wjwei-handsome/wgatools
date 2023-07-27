use crate::parser::chain::ChainReader;
use crate::parser::common::FileFormat;
use crate::parser::maf::MAFReader;
use crate::parser::paf::PAFReader;
use log::{error, warn};
use std::error::Error;
use std::fs::File;
use std::io::{stdin, stdout, BufRead, BufReader, BufWriter, Stdin, Write};
use std::path::Path;

// refer from rustybam, thanks @Mitchell R. Vollger <mrvollger@gmail.com>
type DynResult<T> = Result<T, Box<dyn Error + 'static>>;
const BUFFER_SIZE: usize = 32 * 1024;

/// rational stdin reader: if stdin is empty, exit with error
fn stdin_reader() -> Stdin {
    // check if stdin is empty
    if atty::is(atty::Stream::Stdin) {
        error!("no input content detected");
        std::process::exit(1);
    } else {
        stdin()
    }
}

/// Get a buffer-read input reader from stdin or a file
fn get_input_reader(input: &Option<String>) -> DynResult<Box<dyn BufRead + Send + 'static>> {
    let reader: Box<dyn BufRead + Send + 'static> = match input {
        Some(path) => {
            if path == "-" {
                Box::new(BufReader::with_capacity(BUFFER_SIZE, stdin_reader()))
            } else {
                Box::new(BufReader::with_capacity(BUFFER_SIZE, File::open(path)?))
            }
        }
        None => Box::new(BufReader::with_capacity(BUFFER_SIZE, stdin_reader())),
    };
    Ok(reader)
}

/// get a output writer including stdout and file writer
// TODO: handle IO error
pub fn output_writer(outputpath: &str) -> Box<dyn Write> {
    if outputpath == "-" {
        Box::new(stdout())
    } else {
        let file = File::create(outputpath).unwrap();
        Box::new(BufWriter::new(file))
    }
}

/// check if output file exists and if rewrite
fn outfile_exist(output_file: &String, rewrite: bool) {
    // check if output file exists
    if output_file != "-" {
        let path = Path::new(output_file);
        if path.exists() {
            if rewrite {
                // rewrite the file
                warn!("file {} exist, will rewrite it", output_file);
            } else {
                // exit
                error!("file {} exist, use -r to rewrite it", output_file);
                std::process::exit(1);
            }
        }
    }
}

pub fn convert(
    in_format: FileFormat,
    out_format: FileFormat,
    input: &Option<String>,
    output: &String,
    target_fa_path: &Option<&str>,
    query_fa_path: &Option<&str>,
    rewrite: bool,
) {
    outfile_exist(output, rewrite);
    let reader = get_input_reader(input).expect("Error: cannot read input file");
    match in_format {
        FileFormat::Maf => {
            let mut mafreader = MAFReader::new(reader);
            mafreader.convert(output, out_format);
        }
        FileFormat::Paf => {
            let mut pafreader = PAFReader::new(reader);
            pafreader.convert(output, out_format, *target_fa_path, *query_fa_path);
        }
        FileFormat::Chain => {
            let mut chainreader = ChainReader::new(reader);
            chainreader.convert(output, out_format, *target_fa_path, *query_fa_path);
        }
        _ => {
            error!("Error: unsupported input format");
            std::process::exit(1);
        }
    }
}

/// Command: maf2paf
pub fn maf2paf(input: &Option<String>, output: &String, rewrite: bool) {
    convert(
        FileFormat::Maf,
        FileFormat::Paf,
        input,
        output,
        &None,
        &None,
        rewrite,
    );
}

/// Command: maf2chain
pub fn maf2chain(input: &Option<String>, output: &String, rewrite: bool) {
    convert(
        FileFormat::Maf,
        FileFormat::Chain,
        input,
        output,
        &None,
        &None,
        rewrite,
    );
}

/// Command: paf2chain
pub fn paf2chain(input: &Option<String>, output: &String, rewrite: bool) {
    convert(
        FileFormat::Paf,
        FileFormat::Chain,
        input,
        output,
        &None,
        &None,
        rewrite,
    );
}

/// Command: paf2maf
pub fn paf2maf(
    input: &Option<String>,
    output: &String,
    target_fa_path: &str,
    query_fa_path: &str,
    rewrite: bool,
) {
    let target_fa_path = Some(target_fa_path);
    let query_fa_path = Some(query_fa_path);
    convert(
        FileFormat::Paf,
        FileFormat::Maf,
        input,
        output,
        &target_fa_path,
        &query_fa_path,
        rewrite,
    );
}

/// Command: chain2maf
pub fn chain2maf(
    input: &Option<String>,
    output: &String,
    target_fa_path: &str,
    query_fa_path: &str,
    rewrite: bool,
) {
    let target_fa_path = Some(target_fa_path);
    let query_fa_path = Some(query_fa_path);
    convert(
        FileFormat::Chain,
        FileFormat::Maf,
        input,
        output,
        &target_fa_path,
        &query_fa_path,
        rewrite,
    );
}

/// Command: chain2paf
pub fn chain2paf(input: &Option<String>, output: &String, rewrite: bool) {
    convert(
        FileFormat::Chain,
        FileFormat::Paf,
        input,
        output,
        &None,
        &None,
        rewrite,
    );
}
