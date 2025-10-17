use crate::{
    cli::Cli,
    converter::{chain2maf, chain2paf, maf2chain, maf2paf, maf2sam, paf2chain, paf2maf},
    errors::WGAError,
    parser::{
        chain::ChainReader,
        common::{DotplotMode, DotplotoutFormat, FileFormat},
        maf::MAFReader,
        paf::PAFReader,
    },
    tools::{
        caller::{call_var_maf, call_var_paf},
        chunk::chunk_maf,
        dotplot::dotplot,
        filter::{filter_chain, filter_maf, filter_paf, filter_paf_align_pair},
        index::{build_index, MafIndex},
        mafextra::maf_extract_idx,
        pafcov::pafcov,
        pseudomaf::generate_pesudo_maf,
        rename::rename_maf,
        stat::{stat_maf, stat_paf}, // trimovp::trim_ovp,
        validate::parallel_validatepaf,
    },
};
use clap::CommandFactory;
use clap_complete::{generate, Shell};
use log::{info, warn};
use regex::Regex;
use std::io::{stdin, stdout, BufRead, BufReader, BufWriter, Read, Stdin, Write};
use std::path::Path;
use std::{fs::File, path::PathBuf};

// TODO : define a pub type WResult = Result<(), WGAError>;

const BUFFER_SIZE: usize = 32 * 1024;

const MAGIC_MAX_LEN: usize = 6;
// compressed file magic number, ref: https://docs.rs/infer/latest/infer/archive/index.html
const GZ_MAGIC: [u8; 3] = [0x1f, 0x8b, 0x08];
const BZ_MAGIC: [u8; 3] = [0x42, 0x5a, 0x68];
const XZ_MAGIC: [u8; 6] = [0xfd, 0x37, 0x7a, 0x58, 0x5A, 0x00];

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

fn get_magic_num(path: &str) -> Result<[u8; MAGIC_MAX_LEN], WGAError> {
    let mut buffer: [u8; MAGIC_MAX_LEN] = [0; MAGIC_MAX_LEN];
    let mut fp = File::open(path)?;
    let _ = fp.read(&mut buffer)?;
    Ok(buffer)
}

fn is_gzipped(path: &str) -> Result<bool, WGAError> {
    let buffer = get_magic_num(path)?;
    let gz_or_not =
        buffer[0] == GZ_MAGIC[0] && buffer[1] == GZ_MAGIC[1] && buffer[2] == GZ_MAGIC[2];
    Ok(gz_or_not || Path::new(path).extension().is_some_and(|ext| ext == "gz"))
}

fn is_bzipped(path: &str) -> Result<bool, WGAError> {
    let buffer = get_magic_num(path)?;
    let bz_or_not =
        buffer[0] == BZ_MAGIC[0] && buffer[1] == BZ_MAGIC[1] && buffer[2] == BZ_MAGIC[2];
    Ok(bz_or_not || Path::new(path).extension().is_some_and(|ext| ext == "bz2"))
}

fn is_xz(path: &str) -> Result<bool, WGAError> {
    let buffer = get_magic_num(path)?;
    let xz_or_not = buffer[0] == XZ_MAGIC[0]
        && buffer[1] == XZ_MAGIC[1]
        && buffer[2] == XZ_MAGIC[2]
        && buffer[3] == XZ_MAGIC[3]
        && buffer[4] == XZ_MAGIC[4]
        && buffer[5] == XZ_MAGIC[5];
    Ok(xz_or_not || Path::new(path).extension().is_some_and(|ext| ext == "xz"))
}

pub fn get_input_reader(input: &Option<String>) -> Result<Box<dyn BufRead + Send>, WGAError> {
    let reader: Box<dyn BufRead + Send> = if let Some(path) = input {
        match File::open(path) {
            Ok(file) => {
                if is_xz(path)? {
                    // decode xz compressed file
                    Box::new(BufReader::with_capacity(
                        BUFFER_SIZE,
                        xz2::read::XzDecoder::new_multi_decoder(file),
                    ))
                } else if is_gzipped(path)? {
                    // decode gzip compressed file
                    Box::new(BufReader::with_capacity(
                        BUFFER_SIZE,
                        flate2::read::MultiGzDecoder::new(file),
                    ))
                } else if is_bzipped(path)? {
                    // decode bzip2 compressed file
                    Box::new(BufReader::with_capacity(
                        BUFFER_SIZE,
                        bzip2::read::MultiBzDecoder::new(file),
                    ))
                } else {
                    // stdin flag "-" covered
                    Box::new(BufReader::with_capacity(BUFFER_SIZE, file))
                }
            }
            Err(_) => return Err(WGAError::FileNotExist(PathBuf::from(path))),
        }
    } else {
        Box::new(BufReader::with_capacity(BUFFER_SIZE, stdin_reader()?))
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
    // check if output file exists
    check_outfile(outputpath, rewrite)?;

    // if output is stdout, return stdout writer directly
    if outputpath == "-" {
        return Ok(Box::new(BufWriter::new(stdout())));
    }

    let file = File::create(outputpath)?;
    // TODO: add compression level option
    let compression_level: u32 = 6;

    let writer: Box<dyn Write> = if Path::new(outputpath)
        .extension()
        .is_some_and(|ext| ext == "xz")
    {
        // encode file to xz format
        Box::new(BufWriter::with_capacity(
            BUFFER_SIZE,
            xz2::write::XzEncoder::new(file, compression_level),
        ))
    } else if Path::new(outputpath)
        .extension()
        .is_some_and(|ext| ext == "gz")
    {
        // encode file to gzip format
        Box::new(BufWriter::with_capacity(
            BUFFER_SIZE,
            flate2::write::GzEncoder::new(file, flate2::Compression::new(compression_level)),
        ))
    } else if Path::new(outputpath)
        .extension()
        .is_some_and(|ext| ext == "bz2")
    {
        // encode file to bzip2 format
        Box::new(BufWriter::with_capacity(
            BUFFER_SIZE,
            bzip2::write::BzEncoder::new(file, bzip2::Compression::new(compression_level)),
        ))
    } else if outputpath != "-" {
        Box::new(BufWriter::with_capacity(BUFFER_SIZE, file))
    } else {
        Box::new(BufWriter::new(stdout()))
    };

    Ok(writer)
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
pub fn wrap_maf2paf(
    input: &Option<String>,
    output: &str,
    query_name: Option<String>,
    rewrite: bool,
) -> Result<(), WGAError> {
    // prepare reader and writer
    let (reader, mut writer) = prepare_rdr_wtr(input, output, rewrite)?;
    let mut mafrdr = MAFReader::new(reader)?;
    maf2paf(&mut mafrdr, &mut writer, query_name.as_deref())?;
    Ok(())
}

/// Command: maf2chain
pub fn wrap_maf2chain(
    input: &Option<String>,
    output: &str,
    rewrite: bool,
    query_name: Option<String>,
) -> Result<(), WGAError> {
    // prepare reader and writer
    let (reader, mut writer) = prepare_rdr_wtr(input, output, rewrite)?;
    let mut mafrdr = MAFReader::new(reader)?;
    maf2chain(&mut mafrdr, &mut writer, query_name.as_deref())?;
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
#[allow(clippy::too_many_arguments)]
pub fn wrap_maf_call(
    input: &Option<String>,
    output: &str,
    rewrite: bool,
    snp: bool,
    inv: bool,
    svlen: u64,
    between: bool,
    sample: Option<&str>,
    query_name: Option<&str>,
    query_regex: Option<&Regex>,
    chunk_size: Option<usize>,
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
        inv,
        svlen,
        between,
        sample,
        query_name,
        query_regex,
        chunk_size,
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
/// Command: paf call
pub fn wrap_paf_call(
    input: &Option<String>,
    t_fa_path: &str,
    q_fa_path: &str,
    output: &str,
    rewrite: bool,
    snp: bool,
    svlen: u64,
    between: bool,
    sample: Option<&str>,
) -> Result<(), WGAError> {
    // prepare reader and writer
    let (reader, mut writer) = prepare_rdr_wtr(input, output, rewrite)?;

    // check if fasta files exist
    if !Path::new(t_fa_path).exists() {
        return Err(WGAError::FileNotExist(PathBuf::from(t_fa_path)));
    }
    if !Path::new(q_fa_path).exists() {
        return Err(WGAError::FileNotExist(PathBuf::from(q_fa_path)));
    }

    // check if fasta index files exist, if not create them
    if !Path::new(&format!("{}.fai", t_fa_path)).exists() {
        return Err(WGAError::FileNotExist(PathBuf::from(format!(
            "{}.fai",
            t_fa_path
        ))));
    }
    if !Path::new(&format!("{}.fai", q_fa_path)).exists() {
        return Err(WGAError::FileNotExist(PathBuf::from(format!(
            "{}.fai",
            q_fa_path
        ))));
    }

    // initialize PAF reader
    let mut pafreader = PAFReader::new(reader);

    call_var_paf(
        &mut pafreader,
        t_fa_path,
        q_fa_path,
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
    query_name: Option<String>,
    rewrite: bool,
    each: bool,
) -> Result<(), WGAError> {
    // prepare reader and writer
    let (reader, mut writer) = prepare_rdr_wtr(input, output, rewrite)?;

    // match format and call stat
    match format {
        FileFormat::Maf => {
            let mafrdr = MAFReader::new(reader)?;
            stat_maf(mafrdr, &mut writer, each, query_name.as_deref())?
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
    min_align_size: Option<u64>,
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
            match min_align_size {
                Some(min_align_size) => {
                    warn!("`min_align_size` is set, will not filter paf `min_block_size` and `min_query_size`");
                    filter_paf_align_pair(pafrdr, &mut writer, min_align_size)?
                }
                None => filter_paf(pafrdr, &mut writer, min_block_size, min_query_size)?,
            }
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
    // get input name for INFO
    let input_name = match input {
        Some(path) => path,
        None => "stdin",
    };
    info!("start read file: `{}`", input_name);

    info!("start write file to dir: `{}`", output);
    if output == "-" {
        return Err(WGAError::StdoutNotAllowed);
    }
    // judge if outputdir if exists
    let outputdir = Path::new(output);
    if !outputdir.exists() {
        std::fs::create_dir_all(outputdir)?;
    } else {
        // judge if outputdir is dir
        if !outputdir.is_dir() {
            return Err(WGAError::NotDir(outputdir.to_path_buf()));
        }
        // if rewrite
        if rewrite {
            warn!("output dir `{}` exists, will rewrite it", output);
        } else {
            return Err(WGAError::FileReWrite(output.to_string()));
        }
    }
    // get a reader
    let reader = get_input_reader(input)?;
    let pafrdr = PAFReader::new(reader);
    generate_pesudo_maf(pafrdr, output, fa_path, target)?;
    Ok(())
}

// /// A wrapper for PAF trim overlap
// pub fn wrap_paf_trim_overlap(
//     input: &Option<String>,
//     output: &str,
//     rewrite: bool,
// ) -> Result<(), WGAError> {
//     let (reader, mut writer) = prepare_rdr_wtr(input, output, rewrite)?;
//     let pafrdr = PAFReader::new(reader);
//     trim_ovp(pafrdr, &mut writer)?;
//     Ok(())
// }

/// A wrapper for chunk sub-cmd
pub fn wrap_chunk(
    input: &Option<String>,
    output: &str,
    rewrite: bool,
    length: u64,
) -> Result<(), WGAError> {
    // check length > 0
    if length == 0 {
        return Err(WGAError::Other(anyhow::anyhow!(
            "`length` should be greater than 0"
        )));
    }

    // prepare reader and writer
    let (reader, mut writer) = prepare_rdr_wtr(input, output, rewrite)?;

    let mafrdr = MAFReader::new(reader)?;

    // mafrdr.chunk(&mut writer, chunk_count, chunk_length)?;
    chunk_maf(mafrdr, length, &mut writer)?;
    Ok(())
}

/// A wrapper for dotplot sub-cmd
#[allow(clippy::too_many_arguments)]
pub fn wrap_dotplot(
    input: &Option<String>,
    format: FileFormat,
    out_format: DotplotoutFormat,
    mode: DotplotMode,
    no_identity: bool,
    cutoff: Option<usize>,
    output: &str,
    query_name: Option<String>,
    rewrite: bool,
    color: Option<String>,
) -> Result<(), WGAError> {
    // prepare reader and writer
    let (reader, mut writer) = prepare_rdr_wtr(input, output, rewrite)?;
    // let mafrdr = MAFReader::new(reader)?;
    match mode {
        DotplotMode::BaseLevel => {
            if no_identity {
                warn!("`no_identity` is set, but it's not supported in `BaseLevel` mode");
            }
        }
        DotplotMode::Overview => {
            if cutoff.is_some() {
                warn!("`cutoff` is set, but it's not supported in `Overview` mode");
            }
        }
    }

    // set default cutoff to 50
    let cutoff = cutoff.unwrap_or(50);

    dotplot(
        reader,
        &mut writer,
        format,
        out_format,
        mode,
        no_identity,
        cutoff,
        query_name.as_deref(),
        color.as_deref(),
    )?;
    Ok(())
}

/// A wrapper for gen-auto-completion
pub fn wrap_gencomp(shell: Shell, output: &str, rewrite: bool) -> Result<(), WGAError> {
    let mut cmd = Cli::command();
    let mut writer = get_output_writer(output, rewrite)?;
    generate(shell, &mut cmd, "wgatools", &mut writer);
    Ok(())
}

pub fn wrap_validate(
    input: &Option<String>,
    fix: &Option<String>,
    output: &str,
    rewrite: bool,
) -> Result<(), WGAError> {
    // prepare reader and writer
    let (reader, mut writer) = prepare_rdr_wtr(input, output, rewrite)?;
    let pafrdr = PAFReader::new(reader);

    let fix_writer = match fix {
        Some(path) => {
            warn!("`fix` is set, will try to fix the query|target postion of paf file.It does NOT represent the alignment behavior.");
            if path == "-" {
                warn!("STDOUT mixed the validation information and new fixed paf records");
            }
            let input_path = match input {
                Some(path) => path,
                None => "stdin",
            };
            if path == input_path {
                return Err(WGAError::Other(anyhow::anyhow!(
                    "fixed file should not be the same as output file"
                )));
            }
            let fix_writer = get_output_writer(path, true)?;
            Some(fix_writer)
        }
        None => None,
    };

    let fix_flag = fix.is_some();
    parallel_validatepaf(pafrdr, &mut writer, fix_writer, fix_flag)?;

    Ok(())
}
