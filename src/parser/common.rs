use serde::{Deserialize, Serialize};
use std::fmt;

/// Enum the file types
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum FileFormat {
    Maf,
    Sam,
    Bam,
    Paf,
    Delta,
    Chain,
    Bedpe,
    Unknown,
    Blocks,
}

/// Represented in:
/// - PAF 1-9 columns
/// - CHAIN header lines
/// - MAF header lines
#[derive(Debug, Serialize, Deserialize)]
pub struct SeqInfo {
    pub name: String,
    pub size: u64,
    pub strand: Strand,
    pub start: u64,
    pub end: u64,
}

#[derive(Debug, Clone, PartialEq, Copy, Serialize, Deserialize)]
pub enum Strand {
    #[serde(rename = "+")]
    Positive,
    #[serde(rename = "-")]
    Negative,
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Strand::Positive => write!(f, "+"),
            Strand::Negative => write!(f, "-"),
        }
    }
}

/// Define an alignment block
#[derive(Debug, Copy, Clone, Serialize)]
pub struct Block<'a> {
    pub(crate) query_name: &'a str,
    pub(crate) query_start: u64,
    pub(crate) query_end: u64,
    pub(crate) target_name: &'a str,
    pub(crate) target_start: u64,
    pub(crate) target_end: u64,
    pub(crate) strand: Strand,
}

pub trait AlignRecord {
    fn query_name(&self) -> &str;
    fn query_length(&self) -> u64;
    fn query_start(&self) -> u64;
    fn query_end(&self) -> u64;
    fn query_strand(&self) -> Strand;
    fn target_name(&self) -> &str;
    fn target_length(&self) -> u64;
    fn target_start(&self) -> u64;
    fn target_end(&self) -> u64;
    fn target_strand(&self) -> Strand;
    fn get_cigar_bytes(&self) -> &[u8];
}
