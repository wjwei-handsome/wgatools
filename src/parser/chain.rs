use crate::converter::chain2paf::chain2paf;
use crate::errors::ParseError;
use crate::parser::cigar::parse_chain_to_cigar;
use crate::parser::common::{AlignRecord, FileFormat, SeqInfo, Strand};
use crate::parser::maf::MAFRecord;
use crate::parser::paf::PafRecord;
use nom::bytes::complete::{is_not, tag, take_while};
use nom::character::complete::{line_ending, not_line_ending};
use nom::multi::fold_many1;
use nom::sequence::terminated;
use nom::IResult;
use std::fs::File;
use std::io::{BufReader, Read};
use std::{fmt, io};
use crate::converter::chain2maf::chain2maf;

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
    pub fn records(&mut self) -> ChainRecords {
        let mut data = String::with_capacity(512);
        self.inner.read_to_string(&mut data).unwrap();
        ChainRecords { inner: data }
    }

    /// convert method
    pub fn convert(
        &mut self,
        outputpath: &str,
        format: FileFormat,
        t_fa_path: Option<&str>,
        q_fa_path: Option<&str>,
    ) {
        match format {
            FileFormat::Paf => {
                chain2paf(self, outputpath);
            }
            FileFormat::Maf => {
                chain2maf(self, outputpath, t_fa_path, q_fa_path);
            }
            _ => {}
        }
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
    pub header: Header,
    pub lines: Vec<ChainDataLine>,
}

pub struct ChainRecords {
    inner: String,
}

impl Iterator for ChainRecords {
    type Item = Result<ChainRecord, String>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.inner.is_empty() {
            return None;
        }
        match chain_parser(&self.inner) {
            Ok((i, r)) => {
                self.inner = i.to_string();
                Some(Ok(r))
            }
            Err(e) => {
                let mut msg = format!("{:?}", e);
                msg.push_str(&self.inner);
                Some(Err(msg))
            }
        }
    }
}

/// Define a chain header
#[derive(Debug)]
pub struct Header {
    score: f64, // could be u64?
    target: SeqInfo,
    query: SeqInfo,
    pub(crate) chain_id: usize,
}

/// Define a chain data lines
#[derive(Debug)]
pub struct ChainDataLine {
    pub(crate) size: u64,
    pub(crate) query_diff: u64,
    pub(crate) target_diff: u64,
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

impl From<&PafRecord> for Header {
    fn from(value: &PafRecord) -> Self {
        Header {
            score: value.mapq as f64,
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
        }
    }
}

impl From<&MAFRecord> for Header {
    fn from(value: &MAFRecord) -> Self {
        Header {
            score: 255.0,
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
        }
    }
}

impl fmt::Display for Header {
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

fn parse_header(input: &str) -> Result<Header, ParseError> {
    let mut iter = input.split_whitespace();
    // let _ = iter.next().unwrap(); // skip chain
    let score = iter.next().unwrap().parse::<f64>().unwrap();
    let target_name = iter.next().unwrap();
    let target_size = iter.next().unwrap().parse::<u64>().unwrap();
    let target_strand = iter.next().unwrap().parse::<Strand>().unwrap();
    let target_start = iter.next().unwrap().parse::<u64>().unwrap();
    let target_end = iter.next().unwrap().parse::<u64>().unwrap();
    let query_name = iter.next().unwrap();
    let query_size = iter.next().unwrap().parse::<u64>().unwrap();
    let query_strand = iter.next().unwrap().parse::<Strand>().unwrap();
    let query_start = iter.next().unwrap().parse::<u64>().unwrap();
    let query_end = iter.next().unwrap().parse::<u64>().unwrap();
    let chain_id = iter.next().unwrap().parse::<usize>().unwrap();
    Ok(Header {
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

fn line_not_chain(i: &str) -> IResult<&str, &str> {
    terminated(is_not("chain\n"), line_ending)(i)
}

fn parse_chain_data_line(input: &str) -> IResult<&str, Vec<ChainDataLine>> {
    fold_many1(
        line_not_chain,
        Vec::new,
        |mut acc: Vec<ChainDataLine>, item: &str| {
            let mut dataline = item.split_whitespace();
            let size = dataline.next().unwrap().parse::<u64>().unwrap();
            let query_diff = match dataline.next() {
                Some(query_diff) => query_diff.parse::<u64>().unwrap(),
                None => 0u64,
            };
            let target_diff = match dataline.next() {
                Some(target_diff) => target_diff.parse::<u64>().unwrap(),
                None => 0u64,
            };
            acc.push(ChainDataLine {
                size,
                query_diff,
                target_diff,
            });
            acc
        },
    )(input)
}

fn chain_parser(input: &str) -> IResult<&str, ChainRecord> {
    let (input, _) = tag("chain")(input)?;
    let (input, header_line) = not_line_ending(input)?;
    let header = parse_header(header_line).unwrap();
    let (input, _) = line_ending(input)?;
    let (input, blocks) = parse_chain_data_line(input)?;
    let (input, _) = take_while(|x| x != 'c')(input)?; // should better
    let chainrecord = ChainRecord {
        header,
        lines: blocks,
    };
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

    fn convert2paf(&self) -> PafRecord {
        let cigar = parse_chain_to_cigar(self, false);
        let cigar_string = cigar.cigar_string;
        let block_length =
            (cigar.match_count + cigar.mismatch_count + cigar.ins_count + cigar.del_count) as u64;
        let matches = cigar.match_count as u64;
        PafRecord {
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
        }
    }
}
