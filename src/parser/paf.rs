use crate::errors::WGAError;
use crate::parser::cigar::parse_paf_to_cigar;
use crate::parser::common::{AlignRecord, RecStat, Strand};
use csv::{DeserializeRecordsIter, ReaderBuilder};
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io;
use std::str;

/// Parser for PAF format files
pub struct PAFReader<R: io::Read> {
    inner: csv::Reader<R>,
}

impl<R> PAFReader<R>
where
    R: io::Read + Send,
{
    /// Create a new PAF parser
    pub fn new(reader: R) -> Self {
        PAFReader {
            inner: ReaderBuilder::new()
                .flexible(true)
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

impl PAFReader<File> {
    /// Create a new PAF parser from a file path
    pub fn from_path<P: AsRef<std::path::Path>>(path: P) -> io::Result<PAFReader<File>> {
        File::open(path).map(PAFReader::new)
    }
}

#[derive(Debug, Serialize, Deserialize, Default)]
/// A PAF record refer to https://github.com/lh3/miniasm/blob/master/PAF.md
pub struct PafRecord {
    pub query_name: String,
    pub query_length: u64,
    pub query_start: u64,
    pub query_end: u64,
    pub strand: Strand,
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

/// An iterator struct for PAF records
pub struct Records<'a, R: io::Read> {
    inner: DeserializeRecordsIter<'a, R, PafRecord>,
}

/// impl Iterator for Records
impl<R: io::Read> Iterator for Records<'_, R> {
    type Item = csv::Result<PafRecord>;
    fn next(&mut self) -> Option<csv::Result<PafRecord>> {
        self.inner.next()
    }
}

/// impl AlignRecord Trait for PafRecord
impl AlignRecord for PafRecord {
    fn query_name(&self) -> &str {
        &self.query_name
    }

    fn query_length(&self) -> u64 {
        self.query_length
    }

    fn query_start(&self) -> u64 {
        self.query_start
    }

    fn query_end(&self) -> u64 {
        self.query_end
    }

    fn query_strand(&self) -> Strand {
        self.strand
    }

    fn target_name(&self) -> &str {
        &self.target_name
    }

    fn target_length(&self) -> u64 {
        self.target_length
    }

    fn target_start(&self) -> u64 {
        self.target_start
    }

    fn target_end(&self) -> u64 {
        self.target_end
    }

    fn target_strand(&self) -> Strand {
        Strand::Positive
    }

    fn get_cigar_string(&self) -> Result<String, WGAError> {
        let cg_tag = self.tags.iter().find(|x| x.starts_with("cg:Z:"));
        let cs_tag = self.tags.iter().find(|x| x.starts_with("cs:Z:"));

        match cg_tag {
            Some(cg) => Ok(cg.to_string()),
            None => match cs_tag {
                Some(cs) => {
                    // remove the prefix cs:Z:
                    let cs = &cs[5..];
                    let mut cg_tag = cs_to_cigar(cs);
                    // add cg:Z: prefix
                    cg_tag.insert_str(0, "cg:Z:");
                    Ok(cg_tag)
                }
                None => Err(WGAError::CigarTagNotFound),
            },
        }
    }

    fn target_align_size(&self) -> u64 {
        // self.block_length
        self.target_end - self.target_start
    }

    fn get_stat(&self) -> Result<RecStat, WGAError> {
        // just convert cigar to stat
        let cigar = parse_paf_to_cigar(self)?;
        Ok(RecStat::from(cigar))
    }
}

/// cstag is represented as :6-ata:10+gtc:4*at:3, where :[0-9]+ represents an identical block, -ata represents a deletion, +gtc an insertion and *at indicates reference base a is substituted with a query base t.
/// cgtag : 6M3D10M3I4M1X3M
///     let cs_tag = ":6-ata:10+gtc:4*at*tg:3";
///     let cigar = cs_to_cigar(cs_tag);
///     println!("{}", cigar);  // output: 6M3D10M3I4M2X3M
fn cs_to_cigar(cs_tag: &str) -> String {
    let re = Regex::new(r"(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)").unwrap();
    let mut cigar = String::new();
    let mut last_op = 'M';
    let mut last_len = 0;

    for cap in re.captures_iter(cs_tag) {
        let part = &cap[0];
        match &part[0..1] {
            ":" => {
                let length: usize = part[1..].parse().unwrap();
                if last_op == 'M' {
                    last_len += length;
                } else {
                    if last_len > 0 {
                        cigar.push_str(&format!("{}{}", last_len, last_op));
                    }
                    last_op = 'M';
                    last_len = length;
                }
            }
            "-" => {
                let length = part[1..].len();
                if last_len > 0 {
                    cigar.push_str(&format!("{}{}", last_len, last_op));
                }
                cigar.push_str(&format!("{}D", length));
                last_len = 0;
                last_op = 'M';
            }
            "+" => {
                let length = part[1..].len();
                if last_len > 0 {
                    cigar.push_str(&format!("{}{}", last_len, last_op));
                }
                cigar.push_str(&format!("{}I", length));
                last_len = 0;
                last_op = 'M';
            }
            "*" => {
                if last_op == 'X' {
                    last_len += 1;
                } else {
                    if last_len > 0 {
                        cigar.push_str(&format!("{}{}", last_len, last_op));
                    }
                    last_op = 'X';
                    last_len = 1;
                }
            }
            _ => {}
        }
    }

    if last_len > 0 {
        cigar.push_str(&format!("{}{}", last_len, last_op));
    }

    cigar
}
