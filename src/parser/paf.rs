use csv::{DeserializeRecordsIter, ReaderBuilder};
use nom::bytes::complete::{tag, take_while};
use nom::character::{is_alphabetic, is_digit};
use nom::IResult;
use std::fs::File;
use std::io;
use std::str;

use nom::error::Error;
use nom::multi::fold_many0;
use serde::{Deserialize, Serialize};

/// Parser for PAF format files
pub struct PafReader<R: io::Read> {
    inner: csv::Reader<R>,
}

impl<R> PafReader<R>
where
    R: io::Read,
{
    /// Create a new PAF parser
    pub fn new(reader: R) -> Self {
        PafReader {
            inner: ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .comment(Some(b'#'))
                .from_reader(reader),
        }
    }

    /// Iterate over the records in the PAF file
    pub fn records(&mut self) -> Records<'_, R> {
        Records {
            inner: self.inner.deserialize(),
        }
    }
}

impl PafReader<File> {
    /// Create a new PAF parser from a file path
    pub fn from_path<P: AsRef<std::path::Path>>(path: P) -> io::Result<PafReader<File>> {
        File::open(path).map(PafReader::new)
    }
}

/// An iterator struct for PAF records
pub struct Records<'a, R: io::Read> {
    inner: DeserializeRecordsIter<'a, R, PafRecord>,
}

impl<'a, R: io::Read> Iterator for Records<'a, R> {
    type Item = csv::Result<PafRecord>;

    fn next(&mut self) -> Option<csv::Result<PafRecord>> {
        self.inner.next()
    }
}

#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
/// A PAF record refer to https://github.com/lh3/miniasm/blob/master/PAF.md
pub struct PafRecord {
    pub query_name: String,
    pub query_length: u64,
    pub query_start: u64,
    pub query_end: u64,
    pub strand: char,
    pub target_name: String,
    pub target_length: u64,
    pub target_start: u64,
    pub target_end: u64,
    pub matches: u64,
    pub block_length: u64,
    pub mapq: u64,
    #[serde(default)]
    pub tags: Vec<String>,
}

fn parse_cigar(input: &[u8]) -> IResult<&[u8], (&[u8], &[u8])> {
    if input.is_empty() {
        return Err(nom::Err::Error(Error::new(
            input,
            nom::error::ErrorKind::Eof,
        )));
    }
    println!("{:?}", input);
    let (input, len) = take_while(is_digit)(input)?;
    println!("{:?}", input);
    let (input, op) = take_while(is_alphabetic)(input)?;
    println!("{:?}", input);
    Ok((input, (op, len)))
}

/// "cg:Z:5M1D2I10M" => [("M", 5), ("D", 1), ("I", 2), ("M", 10)]
#[allow(clippy::type_complexity)]
pub fn parse_cigar_to_alignment(cigar: &[u8]) -> IResult<&[u8], Vec<(&[u8], &[u8])>> {
    let (cigar, _tag) = tag(b"cg:Z:")(cigar)?;
    // let mut res_vec = Vec::new();
    let mut qpos = 10u64;
    let mut tpos = 20u64;
    let (rest, res) = fold_many0(
        parse_cigar,
        Vec::new,
        |mut acc: Vec<(&[u8], &[u8])>, item| {
            println!("{:?}", item);
            acc.push(item);
            let count = item.1;
            let op = item.0;
            let count_u64 = str::from_utf8(count).unwrap().parse::<u64>().unwrap();
            println!("count: {}", count_u64);
            match op {
                b"M" => {
                    qpos += count_u64;
                    tpos += count_u64;
                }
                b"I" => {
                    qpos += count_u64;
                }
                b"D" => {
                    tpos += count_u64;
                }
                _ => {}
            };
            println!("qpos: {}, tpos: {}", qpos, tpos);
            acc
        },
    )(cigar)?;
    Ok((rest, res))
}
