use crate::parser::chain::ChainReader;
use crate::parser::common::FileFormat;
use crate::parser::maf::MAFReader;
use crate::parser::paf::PAFReader;
use crate::tools::caller::call_var_maf;
use crate::tools::index::{build_index, MafIndex};
use crate::tools::mafextra::maf_extract_idx;
use log::{error, info, warn};
use std::error::Error;
use std::fs::File;
use std::io::{stdin, stdout, BufRead, BufReader, BufWriter, Stdin, Write};
use std::path::Path;

pub fn reverse_complement(input: &str) -> String {
    let mut output = String::with_capacity(input.len());
    for c in input.chars().rev() {
        match c {
            'A' => output.push('T'),
            'C' => output.push('G'),
            'G' => output.push('C'),
            'T' => output.push('A'),
            'N' => output.push('N'),
            _ => panic!("Invalid character"),
        }
    }
    output
}

// refer from rustybam, thanks @Mitchell R. Vollger <mrvollger@gmail.com>
type DynResult<T> = Result<T, Box<dyn Error + 'static>>;
const BUFFER_SIZE: usize = 32 * 1024;

/// rational stdin reader: if stdin is empty, exit with error
pub fn stdin_reader() -> Stdin {
    // check if stdin is empty
    if atty::is(atty::Stream::Stdin) {
        error!("no input content detected");
        std::process::exit(1);
    } else {
        stdin()
    }
}

/// Get a buffer-read input reader from stdin or a file
pub fn get_input_reader(input: &Option<String>) -> DynResult<Box<dyn BufRead + Send + 'static>> {
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

    let input_name = match input {
        Some(path) => path,
        None => "stdin",
    };

    let output_name = match output.as_str() {
        "-" => "stdout",
        path => path,
    };

    let reader = match get_input_reader(input) {
        Ok(reader) => reader,
        Err(why) => {
            error!("Input Error: {} in {}", why, input_name);
            std::process::exit(1);
        }
    };

    info!("start read {:?} file: {}", in_format, input_name);

    match in_format {
        FileFormat::Maf => {
            let mut mafreader = MAFReader::new(reader);
            info!(
                "start convert {:?} file into {:?}: {}",
                in_format, out_format, output_name
            );
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

/// Command: maf2sam
pub fn maf2sam(input: &Option<String>, output: &String, rewrite: bool) {
    convert(
        FileFormat::Maf,
        FileFormat::Sam,
        input,
        output,
        &None,
        &None,
        rewrite,
    );
}

/// Command: build maf index
pub fn wrap_build_index(input: &String, outputpath: &str) {
    let outputpath = match outputpath {
        "-" => {
            // add .idx suffix to input file
            let mut path = input.clone();
            path.push_str(".index");
            path
        }
        path => path.to_owned(),
    };

    let mut mafreader = MAFReader::from_path(input).unwrap();
    // // init index-writer and csv-writer for deserializing
    // let mut idx_wtr = csv::WriterBuilder::new()
    //     .delimiter(b'\t')
    //     .has_headers(false)
    //     .from_writer(output_writer(&outputpath));
    // init index-writer
    let idx_wtr = output_writer(&outputpath);
    build_index(&mut mafreader, idx_wtr)
}

/// Command: maf extract
pub fn wrap_maf_extract(
    input: &Option<String>,
    regions: &Option<Vec<String>>,
    region_file: &Option<String>,
    output: &String,
    rewrite: bool,
) {
    // judge regions and region_file
    if regions.is_none() && region_file.is_none() {
        error!("regions or region_file must be specified");
        std::process::exit(1);
    }

    outfile_exist(output, rewrite);

    let _input_name = match input {
        Some(path) => path,
        None => "stdin",
    };

    let _output_name = match output.as_str() {
        "-" => "stdout",
        path => path,
    };

    let mut writer = output_writer(output);

    match input {
        // if input if from file, use index
        Some(path) => {
            let mut mafreader = MAFReader::from_path(path).unwrap();
            let index_path = format!("{}.index", path);
            let index_rdr = match File::open(index_path) {
                Ok(file) => BufReader::new(file),
                Err(err) => {
                    warn!(
                        "failed to open index file: {}, please use `maf-index` to create it; will use iterator to extract",
                        err
                    );
                    // std::process::exit(1);
                    todo!("use iterator to extract")
                }
            };
            let mafindex: MafIndex = serde_json::from_reader(index_rdr).unwrap();
            maf_extract_idx(regions, region_file, &mut mafreader, mafindex, &mut writer);
        }
        None => {
            todo!("use iterator to extract")
        }
    }
}

/// Command: maf call
pub fn wrap_maf_call(
    input: &String,
    output: &String,
    rewrite: bool,
    snp: bool,
    svlen: u64,
    between: bool,
    sample: Option<&str>,
) {
    outfile_exist(output, rewrite);

    let _output_name = match output.as_str() {
        "-" => "stdout",
        path => path,
    };

    let mut writer = output_writer(output);

    let mut mafreader = MAFReader::from_path(input).unwrap();
    let index_path = format!("{}.index", input);

    let mafindex = match File::open(index_path) {
        Ok(file) => {
            let index_rdr = BufReader::new(file);
            match serde_json::from_reader(index_rdr) {
                Ok(mafindex) => Some(mafindex),
                Err(err) => {
                    warn!(
                        "failed to deserialize index file: {},
                        please check!",
                        err
                    );
                    None
                }
            }
        }
        Err(err) => {
            warn!(
                "failed to open index file: {},
                please use `maf-index` to create it;
                will not generate contig info",
                err
            );
            None
        }
    };

    call_var_maf(
        &mut mafreader,
        mafindex,
        &mut writer,
        snp,
        svlen,
        between,
        sample,
    );
}
