use crate::converter::paf2block::paf2blocks;
use crate::converter::paf2chain::paf2chain;
use crate::converter::paf2maf::paf2maf;
use crate::errors::WGAError;
use crate::parser::common::{AlignRecord, FileFormat, Strand};
use csv::{DeserializeRecordsIter, ReaderBuilder};
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

    /// convert method
    pub fn convert(
        &mut self,
        outputpath: &str,
        format: FileFormat,
        t_fa_path: Option<&str>,
        q_fa_path: Option<&str>,
    ) {
        match format {
            FileFormat::Chain => paf2chain(self, outputpath),
            FileFormat::Blocks => paf2blocks(self, outputpath),
            FileFormat::Maf => paf2maf(self, outputpath, t_fa_path, q_fa_path),
            _ => {}
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
impl<'a, R: io::Read> Iterator for Records<'a, R> {
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

    fn get_cigar_str(&self) -> Result<&str, WGAError> {
        self.tags
            .iter()
            .find(|x| x.starts_with("cg:Z:"))
            .ok_or(WGAError::CigarTagNotFound)
            .map(|x| x.as_str())
    }

    fn target_align_size(&self) -> u64 {
        self.block_length
    }
}
