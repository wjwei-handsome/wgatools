use crate::parser::maf::MAFRecord;
use crate::parser::paf::PafRecord;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt;
use std::str::FromStr;

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

#[derive(Debug, Clone, PartialEq, Copy, Serialize, Deserialize, Eq)]
pub enum Strand {
    #[serde(rename = "+")]
    Positive,
    #[serde(rename = "-")]
    Negative,
}

impl FromStr for Strand {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Strand::Positive),
            "-" => Ok(Strand::Negative),
            _ => Err(format!("Invalid strand character: {}", s)),
        }
    }
}

impl From<char> for Strand {
    fn from(c: char) -> Self {
        match c {
            '+' => Strand::Positive,
            '-' => Strand::Negative,
            _ => panic!("Invalid strand character"), // TODO: better error handling
        }
    }
}

impl From<&str> for Strand {
    fn from(s: &str) -> Self {
        match s {
            "+" => Strand::Positive,
            "-" => Strand::Negative,
            _ => panic!("Invalid strand character"), // TODO: better error handling
        }
    }
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

/// impl Default for Block
impl<'a> Default for Block<'a> {
    fn default() -> Self {
        Block {
            query_name: "",
            query_start: 0,
            query_end: 0,
            target_name: "",
            target_start: 0,
            target_end: 0,
            strand: Strand::Positive,
        }
    }
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
    fn target_align_size(&self) -> u64;
    fn get_cigar_bytes(&self) -> &[u8] {
        b"*"
    }
    fn get_cigar_string(&self) -> String {
        "*".to_string()
    }
    fn convert2paf(&self) -> PafRecord {
        PafRecord::default()
    }
    fn convert2maf(&self) -> MAFRecord {
        MAFRecord::default()
    }
    fn convert2bam(&self, _name_id_map: &HashMap<&str, u64>) {}
    fn query_seq(&self) -> &str {
        ""
    }
    fn target_seq(&self) -> &str {
        ""
    }
}
