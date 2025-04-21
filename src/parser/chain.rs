use crate::errors::{ParseChainErrKind, WGAError};
use crate::parser::cigar::{parse_chain_to_cigar, parse_cigar_to_trim, parse_maf_seq_to_trim};
use crate::parser::common::{AlignRecord, SeqInfo, Strand};
use crate::parser::maf::MAFRecord;
use crate::parser::paf::PafRecord;
use crate::utils::{parse_str2f64, parse_str2u64};
use nom::bytes::complete::{is_not, tag, take_while};
use nom::character::complete::{line_ending, not_line_ending};
use nom::multi::fold_many1;
use nom::sequence::terminated;
use nom::IResult;
use std::fs::File;
use std::io::{BufReader, Read};
use std::{fmt, io};

/// Reader for MAF file format
pub struct ChainReader<R: Read> {
    inner: BufReader<R>,
}

impl<R> ChainReader<R>
where
    R: Read + Send,
{
    /// Create a new Chain Reader
    pub fn new(reader: R) -> Self {
        ChainReader {
            inner: BufReader::new(reader),
        }
    }

    /// Iterate over the records in the Chain file
    pub fn records(&mut self) -> Result<ChainRecords, WGAError> {
        let mut data = String::with_capacity(512);
        self.inner.read_to_string(&mut data)?;
        Ok(ChainRecords { inner: data })
    }
}

impl ChainReader<File> {
    /// Create a new PAF parser from a file path
    pub fn from_path<P: AsRef<std::path::Path>>(path: P) -> io::Result<ChainReader<File>> {
        File::open(path).map(ChainReader::new)
    }
}

/// Define a chain record and its header and data lines
/// refer into https://genome.ucsc.edu/goldenPath/help/chain.html
#[derive(Debug)]
pub struct ChainRecord {
    pub header: ChainHeader,
    pub lines: Vec<ChainDataLine>,
}

pub struct ChainRecords {
    inner: String,
}

impl Iterator for ChainRecords {
    type Item = Result<ChainRecord, WGAError>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.inner.is_empty() {
            return None;
        }
        match chain_parser(&self.inner) {
            Ok((i, r)) => {
                self.inner = i.to_string();
                Some(Ok(r))
            }
            Err(e) => Some(Err(e)),
        }
    }
}

/// Define a chain header
#[derive(Debug, Default)]
pub struct ChainHeader {
    score: f64, // could be u64?
    target: SeqInfo,
    query: SeqInfo,
    pub chain_id: usize,
}

/// Define a chain data lines
#[derive(Debug, Default)]
pub struct ChainDataLine {
    pub size: u64,
    pub query_diff: u64,
    pub target_diff: u64,
}

impl fmt::Display for ChainDataLine {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\n{}\t{}\t{}",
            self.size, self.query_diff, self.target_diff
        )
    }
}

// we can't use T, refto: https://github.com/rust-lang/rust/issues/50238
impl TryFrom<&MAFRecord> for ChainHeader {
    type Error = WGAError;

    fn try_from(value: &MAFRecord) -> Result<Self, Self::Error> {
        let mut header = ChainHeader {
            score: 255_f64,
            target: SeqInfo {
                name: value.target_name().to_owned(),
                size: value.target_length(),
                strand: Strand::Positive,
                start: value.target_start(),
                end: value.target_end(),
            },
            query: SeqInfo {
                name: value.query_name().to_owned(),
                size: value.query_length(),
                strand: value.query_strand(),
                start: value.query_start(),
                end: value.query_end(),
            },
            chain_id: 0,
        };
        let (head_ins, head_del, tail_ins, tail_del) = parse_maf_seq_to_trim(value)?;
        match value.query_strand() {
            Strand::Positive => {
                header.query.start += head_ins;
                header.target.start += head_del;
                header.query.end -= tail_ins;
                header.target.end -= tail_del;
            }
            Strand::Negative => {
                header.target.start += head_del;
                header.target.end -= tail_del;
                header.query.start = header.query.size - (header.query.end - head_ins);
                header.query.end = header.query.size - (header.query.start + tail_ins);
            }
        }
        Ok(header)
    }
}

impl TryFrom<&PafRecord> for ChainHeader {
    type Error = WGAError;

    fn try_from(value: &PafRecord) -> Result<Self, Self::Error> {
        let mut header = ChainHeader {
            score: 255_f64,
            target: SeqInfo {
                name: value.target_name.clone(),
                size: value.target_length,
                strand: Strand::Positive,
                start: value.target_start,
                end: value.target_end,
            },
            query: SeqInfo {
                name: value.query_name.clone(),
                size: value.query_length,
                strand: value.query_strand(),
                start: value.query_start,
                end: value.query_end,
            },
            chain_id: 0,
        };
        let (head_ins, head_del, tail_ins, tail_del) = parse_cigar_to_trim(value)?;
        match value.query_strand() {
            Strand::Positive => {
                header.query.start += head_ins;
                header.target.start += head_del;
                header.query.end -= tail_ins;
                header.target.end -= tail_del;
            }
            Strand::Negative => {
                header.target.start += head_del;
                header.target.end -= tail_del;
                header.query.start = header.query.size - (header.query.end - head_ins);
                header.query.end = header.query.size - (header.query.start + tail_ins);
            }
        }
        Ok(header)
    }
}

impl fmt::Display for ChainHeader {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "chain\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.score,
            self.target.name,
            self.target.size,
            self.target.strand,
            self.target.start,
            self.target.end,
            self.query.name,
            self.query.size,
            self.query.strand,
            self.query.start,
            self.query.end,
            self.chain_id,
        )
    }
}

fn parse_header(input: &str) -> Result<ChainHeader, WGAError> {
    let mut iter = input.split_whitespace();
    let score = match iter.next() {
        Some(score) => parse_str2f64(score)?,
        None => {
            return Err(WGAError::ParseChain(ParseChainErrKind::FiledMissing(
                "score".to_string(),
            )))
        }
    };
    let target_name = match iter.next() {
        Some(target_name) => target_name,
        None => {
            return Err(WGAError::ParseChain(ParseChainErrKind::FiledMissing(
                "target_name".to_string(),
            )))
        }
    };
    let target_size = match iter.next() {
        Some(target_size) => parse_str2u64(target_size)?,
        None => {
            return Err(WGAError::ParseChain(ParseChainErrKind::FiledMissing(
                "target_size".to_string(),
            )))
        }
    };
    let target_strand = match iter.next() {
        Some(target_strand) => target_strand.parse::<Strand>()?,
        None => {
            return Err(WGAError::ParseChain(ParseChainErrKind::FiledMissing(
                "target_strand".to_string(),
            )))
        }
    };
    let target_start = match iter.next() {
        Some(target_start) => parse_str2u64(target_start)?,
        None => {
            return Err(WGAError::ParseChain(ParseChainErrKind::FiledMissing(
                "target_start".to_string(),
            )))
        }
    };
    let target_end = match iter.next() {
        Some(target_end) => parse_str2u64(target_end)?,
        None => {
            return Err(WGAError::ParseChain(ParseChainErrKind::FiledMissing(
                "target_end".to_string(),
            )))
        }
    };
    let query_name = match iter.next() {
        Some(query_name) => query_name,
        None => {
            return Err(WGAError::ParseChain(ParseChainErrKind::FiledMissing(
                "query_name".to_string(),
            )))
        }
    };
    let query_size = match iter.next() {
        Some(query_size) => parse_str2u64(query_size)?,
        None => {
            return Err(WGAError::ParseChain(ParseChainErrKind::FiledMissing(
                "query_size".to_string(),
            )))
        }
    };
    let query_strand = match iter.next() {
        Some(query_strand) => query_strand.parse::<Strand>()?,
        None => {
            return Err(WGAError::ParseChain(ParseChainErrKind::FiledMissing(
                "query_strand".to_string(),
            )))
        }
    };
    let query_start = match iter.next() {
        Some(query_start) => parse_str2u64(query_start)?,
        None => {
            return Err(WGAError::ParseChain(ParseChainErrKind::FiledMissing(
                "query_start".to_string(),
            )))
        }
    };
    let query_end = match iter.next() {
        Some(query_end) => parse_str2u64(query_end)?,
        None => {
            return Err(WGAError::ParseChain(ParseChainErrKind::FiledMissing(
                "query_end".to_string(),
            )))
        }
    };
    let chain_id = match iter.next() {
        Some(chain_id) => parse_str2u64(chain_id)? as usize,
        None => {
            return Err(WGAError::ParseChain(ParseChainErrKind::FiledMissing(
                "chain_id".to_string(),
            )))
        }
    };
    Ok(ChainHeader {
        score,
        target: SeqInfo {
            name: target_name.to_owned(),
            size: target_size,
            strand: target_strand,
            start: target_start,
            end: target_end,
        },
        query: SeqInfo {
            name: query_name.to_owned(),
            size: query_size,
            strand: query_strand,
            start: query_start,
            end: query_end,
        },
        chain_id,
    })
}

// helper function for chain_parser, terminated if meet another "chain"
fn line_not_chain(i: &str) -> IResult<&str, &str> {
    terminated(is_not("chain\n"), line_ending)(i)
}

// parse line to ChainDataLine
fn parse_line_to_cdl(line: &str) -> Result<ChainDataLine, WGAError> {
    let mut dataline = line.split_whitespace();
    let size = parse_str2u64(dataline.next().ok_or(WGAError::ParseChain(
        ParseChainErrKind::FiledMissing("size".to_string()),
    ))?)?;
    let query_diff = match dataline.next() {
        Some(query_diff) => parse_str2u64(query_diff)?,
        None => 0u64,
    };
    let target_diff = match dataline.next() {
        Some(target_diff) => parse_str2u64(target_diff)?,
        None => 0u64,
    };
    Ok(ChainDataLine {
        size,
        query_diff,
        target_diff,
    })
}

// helper function for chain_parser, fold_many1
fn init_res_vec_cdl() -> Result<Vec<ChainDataLine>, WGAError> {
    let res_vec = Vec::new();
    Ok(res_vec)
}

fn parse_chain_data_line(
    input: &str,
) -> Result<(&str, Result<Vec<ChainDataLine>, WGAError>), WGAError> {
    fold_many1(
        line_not_chain,
        init_res_vec_cdl,
        |mut acc: Result<Vec<ChainDataLine>, WGAError>, item| {
            if let Ok(v) = acc.as_mut() {
                v.push(parse_line_to_cdl(item)?)
            }
            acc
        },
    )(input)
    .map_err(|e| e.into())
}

// parse chain record, return rest input and a ChainRecord
fn chain_parser(input: &str) -> Result<(&str, ChainRecord), WGAError> {
    let (input, _) = tag("chain")(input)?;
    let (input, header_line) = not_line_ending(input)?;
    let header = parse_header(header_line)?;
    let (input, _) = line_ending(input)?;
    let (input, lines) = parse_chain_data_line(input)?;
    let lines = lines?;
    let (input, _) = take_while(|x| x != 'c')(input)?; // should better
    let chainrecord = ChainRecord { header, lines };
    Ok((input, chainrecord))
}

impl AlignRecord for ChainRecord {
    fn query_name(&self) -> &str {
        &self.header.query.name
    }

    fn query_length(&self) -> u64 {
        self.header.query.size
    }

    fn query_start(&self) -> u64 {
        self.header.query.start
    }

    fn query_end(&self) -> u64 {
        self.header.query.end
    }

    fn query_strand(&self) -> Strand {
        self.header.query.strand
    }

    fn target_name(&self) -> &str {
        &self.header.target.name
    }

    fn target_length(&self) -> u64 {
        self.header.target.size
    }

    fn target_start(&self) -> u64 {
        self.header.target.start
    }

    fn target_end(&self) -> u64 {
        self.header.target.end
    }

    fn target_strand(&self) -> Strand {
        self.header.target.strand
    }

    fn target_align_size(&self) -> u64 {
        self.header.target.end - self.header.target.start
    }

    fn convert2paf(&mut self, _query_name: Option<&str>) -> Result<PafRecord, WGAError> {
        let cigar = parse_chain_to_cigar(self, false);
        let cigar_string = String::from("cg:Z:") + &cigar.cigar_string;
        let block_length =
            (cigar.match_count + cigar.mismatch_count + cigar.del_count + cigar.inv_del_count)
                as u64;
        let matches = cigar.match_count as u64;
        Ok(PafRecord {
            query_name: self.query_name().to_string(),
            query_length: self.query_length(),
            query_start: self.query_start(),
            query_end: self.query_end(),
            strand: self.query_strand(),
            target_name: self.target_name().to_string(),
            target_length: self.target_length(),
            target_start: self.target_start(),
            target_end: self.target_end(),
            matches,
            block_length,
            mapq: 255,
            tags: vec![cigar_string],
        })
    }
}
