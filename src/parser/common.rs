use crate::parser::maf::MAFRecord;
use crate::parser::paf::PafRecord;
use clap::ValueEnum;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt;
use std::str::FromStr;

use super::cigar::Cigar;

/// Enum the file types
#[derive(Debug, PartialEq, Clone, Copy, ValueEnum)]
pub enum FileFormat {
    Maf,
    Sam,
    // Bam,
    Paf,
    Chain,
    // Bedpe,
    // Unknown,
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

#[derive(Debug, Default)]
pub struct RecStat {
    pub aligned_size: usize,
    pub matched: usize,
    pub mismatched: usize,
    pub ins_event: usize,
    pub del_event: usize,
    pub ins_size: usize,
    pub del_size: usize,
    pub inv_ins_event: usize,
    pub inv_ins_size: usize,
    pub inv_del_event: usize,
    pub inv_del_size: usize,
    pub inv_event: usize,
    pub inv_size: f32,
}

// Statistic for each record by CIGAR
impl From<Cigar> for RecStat {
    fn from(cigar: Cigar) -> Self {
        let mut rec_stat = RecStat::default();
        rec_stat.matched = cigar.match_count;
        rec_stat.mismatched = cigar.mismatch_count;
        rec_stat.ins_event = cigar.ins_event;
        rec_stat.del_event = cigar.del_event;
        rec_stat.ins_size = cigar.ins_count;
        rec_stat.del_size = cigar.del_count;
        rec_stat.inv_ins_event = cigar.inv_ins_event;
        rec_stat.inv_ins_size = cigar.inv_ins_count;
        rec_stat.inv_del_event = cigar.inv_del_event;
        rec_stat.inv_del_size = cigar.inv_del_count;
        rec_stat.aligned_size =
            rec_stat.matched + rec_stat.mismatched + rec_stat.del_size + rec_stat.inv_del_size;
        let query_align_size =
            rec_stat.matched + rec_stat.mismatched + rec_stat.ins_size + rec_stat.inv_ins_size;
        rec_stat.inv_event = cigar.inv_event;
        if rec_stat.inv_event != 0 {
            rec_stat.inv_size =
                (rec_stat.aligned_size + query_align_size) as f32 / (rec_stat.inv_event + 1) as f32;
        };
        rec_stat
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
    fn get_stat(&self) -> RecStat {
        RecStat::default()
    }
}
