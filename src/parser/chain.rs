use crate::parser::common::{AlignRecord, SeqInfo, Strand};
use crate::parser::paf::PafRecord;
use std::fmt;

/// Define a chain record and its header and data lines
/// refer into https://genome.ucsc.edu/goldenPath/help/chain.html
#[derive(Debug)]
pub struct ChainRecord {
    pub header: Header,
    pub lines: Vec<ChainDataLine>,
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
