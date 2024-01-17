use crate::{
    converter::{chain2maf, chain2paf, maf2chain, maf2paf, maf2sam, paf2chain, paf2maf},
    errors::WGAError,
    parser::{chain::ChainReader, common::FileFormat, maf::MAFReader, paf::PAFReader},
    tools::{
        caller::call_var_maf,
        filter::{filter_chain, filter_maf, filter_paf},
        index::{build_index, MafIndex},
        mafextra::maf_extract_idx,
        pafcov::pafcov,
        pseudomaf::generate_pesudo_maf,
        rename::rename_maf,
        stat::{stat_maf, stat_paf},
    },
};
use log::{info, warn};
use std::io::{stdin, stdout, BufRead, BufReader, BufWriter, Stdin, Write};
use std::path::Path;
use std::{fs::File, path::PathBuf};

// TODO : define a pub type WResult = Result<(), WGAError>;

const BUFFER_SIZE: usize = 32 * 1024;

type RdrWtr = (Box<dyn BufRead + Send>, Box<dyn Write>);
fn prepare_rdr_wtr(
    input: &Option<String>,
    output: &str,
    rewrite: bool,
) -> Result<RdrWtr, WGAError> {
    // get input name for INFO
    let input_name = match input {
        Some(path) => path,
        None => "stdin",
    };
    info!("start read file: `{}`", input_name);

    // init writer and check if output file exists
    let writer = get_output_writer(output, rewrite)?;
    let output_name = match output {
        "-" => "stdout",
        path => path,
    };
    info!("start write file: `{}`", output_name);

    // get a reader
    let reader = get_input_reader(input)?;
    Ok((reader, writer))
}

pub fn parse_str2u64(s: &str) -> Result<u64, WGAError> {
    match s.parse::<u64>() {
        Ok(n) => Ok(n),
        Err(_) => Err(WGAError::ParseIntError(s.to_string())),
    }
}

pub fn parse_str2f64(s: &str) -> Result<f64, WGAError> {
    match s.parse::<f64>() {
        Ok(n) => Ok(n),
        Err(_) => Err(WGAError::ParseFloatError(s.to_string())),
    }
}

pub fn reverse_complement(input: &str) -> Result<String, WGAError> {
    let mut output = String::with_capacity(input.len());
    for c in input.chars().rev() {
        match c {
            'A' => output.push('T'),
            'C' => output.push('G'),
            'G' => output.push('C'),
            'T' => output.push('A'),
            'N' => output.push('N'),
            'a' => output.push('t'),
            'c' => output.push('g'),
            'g' => output.push('c'),
            't' => output.push('a'),
            'n' => output.push('n'),
            _ => return Err(WGAError::InvalidBase(c.to_string())),
        }
    }
    Ok(output)
}

pub fn get_input_reader(input: &Option<String>) -> Result<Box<dyn BufRead + Send>, WGAError> {
    let reader: Box<dyn BufRead + Send> = match input {
        Some(path) => {
            if path == "-" {
                Box::new(BufReader::with_capacity(BUFFER_SIZE, stdin_reader()?))
            } else {
                match File::open(path) {
                    Ok(file) => Box::new(BufReader::with_capacity(BUFFER_SIZE, file)),
                    Err(_) => return Err(WGAError::FileNotExist(PathBuf::from(path))),
                }
            }
        }
        None => Box::new(BufReader::with_capacity(BUFFER_SIZE, stdin_reader()?)),
    };
    Ok(reader)
}

/// rational stdin reader: if stdin is empty, exit with error
fn stdin_reader() -> Result<Stdin, WGAError> {
    // check if stdin is empty
    if atty::is(atty::Stream::Stdin) {
        Err(WGAError::EmptyStdin)
    } else {
        Ok(stdin())
    }
}

fn get_output_writer(outputpath: &str, rewrite: bool) -> Result<Box<dyn Write>, WGAError> {
    check_outfile(outputpath, rewrite)?;
    if outputpath == "-" {
        Ok(Box::new(stdout()))
    } else {
        let file = File::create(outputpath)?;
        Ok(Box::new(BufWriter::new(file)))
    }
}

/// check if output file exists and if rewrite
fn check_outfile(output_file: &str, rewrite: bool) -> Result<(), WGAError> {
    // check if output file exists
    if output_file != "-" {
        let path = Path::new(output_file);
        if path.exists() {
            if rewrite {
                // rewrite the file
                warn!("file {} exist, will rewrite it", output_file);
            } else {
                // exit
                return Err(WGAError::FileReWrite(output_file.to_string()));
            }
        }
    }
    Ok(())
}

/// Command: maf2paf
pub fn wrap_maf2paf(input: &Option<String>, output: &str, rewrite: bool) -> Result<(), WGAError> {
    // prepare reader and writer
    let (reader, mut writer) = prepare_rdr_wtr(input, output, rewrite)?;
    let mut mafrdr = MAFReader::new(reader)?;
    maf2paf(&mut mafrdr, &mut writer)?;
    Ok(())
}

/// Command: maf2chain
pub fn wrap_maf2chain(input: &Option<String>, output: &str, rewrite: bool) -> Result<(), WGAError> {
    // prepare reader and writer
    let (reader, mut writer) = prepare_rdr_wtr(input, output, rewrite)?;
    let mut mafrdr = MAFReader::new(reader)?;
    maf2chain(&mut mafrdr, &mut writer)?;
    Ok(())
}

/// Command: maf2sam
pub fn wrap_maf2sam(input: &Option<String>, output: &str, rewrite: bool) -> Result<(), WGAError> {
    // prepare reader and writer
    let (reader, mut writer) = prepare_rdr_wtr(input, output, rewrite)?;
    let mut mafrdr = MAFReader::new(reader)?;
    maf2sam(&mut mafrdr, &mut writer)?;
    Ok(())
}

/// Command: paf2chain
pub fn wrap_paf2chain(input: &Option<String>, output: &str, rewrite: bool) -> Result<(), WGAError> {
    // prepare reader and writer
    let (reader, mut writer) = prepare_rdr_wtr(input, output, rewrite)?;
    let mut pafrdr = PAFReader::new(reader);
    paf2chain(&mut pafrdr, &mut writer)?;
    Ok(())
}

/// Command: paf2maf
pub fn wrap_paf2maf(
    input: &Option<String>,
    output: &str,
    target_fa_path: &str,
    query_fa_path: &str,
    rewrite: bool,
) -> Result<(), WGAError> {
    // prepare reader and writer
    let (reader, mut writer) = prepare_rdr_wtr(input, output, rewrite)?;
    let mut pafrdr = PAFReader::new(reader);
    paf2maf(&mut pafrdr, &mut writer, target_fa_path, query_fa_path)?;
    Ok(())
}

/// Command: chain2maf
pub fn wrap_chain2maf(
    input: &Option<String>,
    output: &str,
    target_fa_path: &str,
    query_fa_path: &str,
    rewrite: bool,
) -> Result<(), WGAError> {
    // prepare reader and writer
    let (reader, mut writer) = prepare_rdr_wtr(input, output, rewrite)?;
    let mut chainrdr = ChainReader::new(reader);
    chain2maf(&mut chainrdr, &mut writer, target_fa_path, query_fa_path)?;
    Ok(())
}

/// Command: chain2paf
pub fn wrap_chain2paf(input: &Option<String>, output: &str, rewrite: bool) -> Result<(), WGAError> {
    // prepare reader and writer
    let (reader, mut writer) = prepare_rdr_wtr(input, output, rewrite)?;
    let mut chainrdr = ChainReader::new(reader);
    chain2paf(&mut chainrdr, &mut writer)?;
    Ok(())
}

/// Command: build maf index
pub fn wrap_build_index(input: &String, outputpath: &str) -> Result<(), WGAError> {
    let outputpath = match outputpath {
        "-" => {
            // add .idx suffix to input file
            let mut path = input.clone();
            path.push_str(".index");
            path
        }
        path => path.to_owned(),
    };

    let mut mafreader = MAFReader::from_path(input)?;

    // NOTE: new index file will always overwrite old one
    let idx_wtr = get_output_writer(&outputpath, true)?;
    build_index(&mut mafreader, idx_wtr)
}

/// Command: maf extract
pub fn wrap_maf_extract(
    input: &Option<String>,
    regions: &Option<Vec<String>>,
    region_file: &Option<String>,
    output: &str,
    rewrite: bool,
) -> Result<(), WGAError> {
    // judge regions and region_file
    if regions.is_none() && region_file.is_none() {
        return Err(WGAError::EmptyRegion);
    }

    // init writer and check if output file exists
    let output_name = match output {
        "-" => "stdout",
        path => path,
    };
    info!("start write file: `{}`", output_name);
    let mut writer = get_output_writer(output, rewrite)?;

    match input {
        // if input if from file, use index
        Some(path) => {
            if path == "-" {
                return Err(WGAError::StdinNotAllowed);
            }
            let mut mafreader = MAFReader::from_path(path)?;
            let index_path = format!("{}.index", path);
            let index_rdr = BufReader::new(File::open(index_path)?);
            let mafindex: MafIndex = serde_json::from_reader(index_rdr)?;
            let failed_regions =
                maf_extract_idx(regions, region_file, &mut mafreader, mafindex, &mut writer)?;
            for region in failed_regions {
                let err = WGAError::FailedRegion(region);
                warn!("{}", err);
            }
            Ok(())
        }
        // if input is from stdin, raise error
        None => Err(WGAError::StdinNotAllowed),
    }
}

/// Command: maf call
pub fn wrap_maf_call(
    input: &Option<String>,
    output: &str,
    rewrite: bool,
    snp: bool,
    svlen: u64,
    between: bool,
    sample: Option<&str>,
) -> Result<(), WGAError> {
    // prepare reader and writer
    let (reader, mut writer) = prepare_rdr_wtr(input, output, rewrite)?;

    // get mafindex if input is not stdin
    let mafindex = match input {
        Some(path) => {
            if path == "-" {
                None
            } else {
                let index_path = format!("{}.index", path);
                match File::open(index_path) {
                    Err(_) => None,
                    Ok(index_file) => {
                        let index_rdr = BufReader::new(index_file);
                        let mafindex: MafIndex = serde_json::from_reader(index_rdr)?;
                        Some(mafindex)
                    }
                }
            }
        }
        None => None,
    };
    if mafindex.is_none() {
        warn!("maf index not found, will not generate contig info");
    }

    // get mafreader
    let mut mafreader = MAFReader::new(reader)?;

    call_var_maf(
        &mut mafreader,
        mafindex,
        &mut writer,
        snp,
        svlen,
        between,
        sample,
    )?;
    Ok(())
}

/// A wrapper for stat sub-cmd, match format and call `stat_{maf,paf}`
pub fn wrap_stat(
    format: FileFormat,
    input: &Option<String>,
    output: &str,
    rewrite: bool,
    each: bool,
) -> Result<(), WGAError> {
    // prepare reader and writer
    let (reader, mut writer) = prepare_rdr_wtr(input, output, rewrite)?;

    // match format and call stat
    match format {
        FileFormat::Maf => {
            let mafrdr = MAFReader::new(reader)?;
            stat_maf(mafrdr, &mut writer, each)?
        }
        FileFormat::Paf => {
            let pafrdr = PAFReader::new(reader);
            stat_paf(pafrdr, &mut writer, each)?
        }
        _ => {
            return Err(WGAError::NotImplemented);
        }
    }
    Ok(())
}

/// A wrapper for filter sub-cmd, match format and call `filter_{maf,paf}`
pub fn wrap_filter(
    format: FileFormat,
    input: &Option<String>,
    output: &str,
    rewrite: bool,
    min_block_size: u64,
    min_query_size: u64,
) -> Result<(), WGAError> {
    // prepare reader and writer
    let (reader, mut writer) = prepare_rdr_wtr(input, output, rewrite)?;

    match format {
        FileFormat::Maf => {
            let mafrdr = MAFReader::new(reader)?;
            filter_maf(mafrdr, &mut writer, min_block_size, min_query_size)?
        }
        FileFormat::Paf => {
            let pafrdr = PAFReader::new(reader);
            filter_paf(pafrdr, &mut writer, min_block_size, min_query_size)?
        }
        FileFormat::Chain => {
            let chainrdr = ChainReader::new(reader);
            filter_chain(chainrdr, &mut writer, min_block_size, min_query_size)?
        }
        _ => {
            return Err(WGAError::NotImplemented);
        }
    }
    Ok(())
}

/// A wrapper for filter sub-cmd, match format and call `filter_{maf,paf}`
pub fn wrap_rename_maf(
    input: &Option<String>,
    output: &str,
    rewrite: bool,
    prefixs: &[String],
) -> Result<(), WGAError> {
    // prepare reader and writer
    let (reader, mut writer) = prepare_rdr_wtr(input, output, rewrite)?;
    let mafrdr = MAFReader::new(reader)?;
    let prefixs = prefixs.iter().map(|s| s.as_str()).collect::<Vec<&str>>();
    rename_maf(mafrdr, &mut writer, prefixs)?;
    Ok(())
}

/// A wrapper for PAF Converage count
pub fn wrap_paf_cov(input: &Option<String>, output: &str, rewrite: bool) -> Result<(), WGAError> {
    let (reader, mut writer) = prepare_rdr_wtr(input, output, rewrite)?;
    let pafrdr = PAFReader::new(reader);
    pafcov(pafrdr, &mut writer)?;
    Ok(())
}

/// A wrapper for PAF pesudo maf
pub fn wrap_paf_pesudo_maf(
    input: &Option<String>,
    output: &str,
    rewrite: bool,
    fa_path: &Option<String>,
    target: &Option<String>,
) -> Result<(), WGAError> {
    let (reader, mut writer) = prepare_rdr_wtr(input, output, rewrite)?;
    let pafrdr = PAFReader::new(reader);
    generate_pesudo_maf(pafrdr, &mut writer, fa_path, target)?;
    Ok(())
}
