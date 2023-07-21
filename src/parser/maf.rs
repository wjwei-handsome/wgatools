use crate::converter::maf2bam::maf2bam;
use crate::converter::maf2chain::maf2chain;
use crate::converter::maf2paf::maf2paf;
use crate::errors::ParseError;
use crate::parser::cigar::parse_maf_seq_to_cigar;
use crate::parser::common::{AlignRecord, FileFormat, Strand};
use crate::parser::paf::PafRecord;
use noodles_core::Position;
use noodles_sam::alignment::Record as SamRecord;
use noodles_sam::record::{Cigar, Data, Flags, ReadName, Sequence};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, Read};

/// Parser for MAF file format
pub struct MAFReader<R: Read> {
    inner: BufReader<R>,
    pub header: String,
}

impl<R> MAFReader<R>
where
    R: Read+Send,
{
    /// Create a new PAF parser
    pub fn new(reader: R) -> Self {
        let mut buf_reader = BufReader::new(reader);
        let mut header = String::new();
        buf_reader.read_line(&mut header).unwrap();
        MAFReader {
            inner: buf_reader,
            header,
        }
    }

    /// Iterate over the records in the MAF file
    /// ```
    ///use wgalib::parser::maf::MAFReader;
    /// let mut reader = MAFReader::from_path("data/test.maf").unwrap();
    ///for record in reader.records() {
    ///    let record = record.unwrap();
    ///    println!("{:?}", record);
    ///}
    /// ```
    pub fn records(&mut self) -> MAFRecords<R> {
        MAFRecords {
            inner: self.inner.by_ref(),
        }
    }

    /// convert method
    pub fn convert(&mut self, outputpath: &str, format: FileFormat) {
        match format {
            FileFormat::Chain => {
                maf2chain(self, outputpath);
            }
            FileFormat::Bam => {
                maf2bam(self, outputpath);
            }
            FileFormat::Paf => {
                maf2paf(self, outputpath);
            }
            _ => {}
        }
    }
}

impl MAFReader<File> {
    /// Create a new PAF parser from a file path
    pub fn from_path<P: AsRef<std::path::Path>>(path: P) -> io::Result<MAFReader<File>> {
        File::open(path).map(MAFReader::new)
    }
}

/// A MAF s-line refer to https://genome.ucsc.edu/FAQ/FAQformat.html#format5
// a score=111
// s ref    100 10 + 100000 ---AGC-CAT-CATT
// s contig 0   10 + 10     ---AGC-CAT-CATT
//
// a score=222
// s ref    100 12 + 100000 ---AGC-CAT-CATTTT
// s contig 0   12 + 12     ---AGC-CAT-CATTTT
#[derive(Debug, PartialEq, Eq)]
struct MAFSLine {
    mode: char,
    name: String,
    start: u64,
    align_size: u64,
    strand: Strand,
    size: u64,
    seq: String,
}

fn str2u64(s: &str) -> Result<u64, ParseError> {
    // TODO: move to common.rs module
    match s.parse::<u64>() {
        Ok(n) => Ok(n),
        Err(_) => Err(ParseError::new_parse_int_err(s)),
    }
}

fn parse_sline(line: String) -> Result<MAFSLine, ParseError> {
    let mut iter = line.split_whitespace();
    let mode = match iter.next() {
        Some(mode) => mode.chars().next().unwrap(), // TODO: error handling
        None => panic!("s-line mode is missing"),   // TODO: error handling
    };
    let name = match iter.next() {
        Some(name) => name.to_string(),
        None => panic!("s-line name is missing"), // TODO: error handling
    };
    let start = match iter.next() {
        Some(start) => str2u64(start)?,
        None => panic!("s-line start is missing"), // TODO: error handling
    };
    let align_size = match iter.next() {
        Some(align_size) => str2u64(align_size)?, // TODO: error handling
        None => panic!("s-line align_size is missing"), // TODO: error handling
    };
    let strand = match iter.next() {
        Some(strand) => Strand::from(strand), // TODO: error handling
        None => panic!("s-line strand is missing"), // TODO: error handling
    };
    let size = match iter.next() {
        Some(size) => str2u64(size)?,
        None => panic!("s-line size is missing"), // TODO: error handling
    };
    let seq = match iter.next() {
        Some(seq) => seq.to_string(),
        None => panic!("s-line seq is missing"), // TODO: error handling
    };
    if iter.next().is_some() {
        panic!("s-line has more than 8 fields")
    };
    Ok(MAFSLine {
        mode,
        name,
        start,
        align_size,
        strand,
        size,
        seq,
    })
}

fn sline_from_string(value: String) -> Result<MAFSLine, ParseError> {
    let s_line = parse_sline(value)?;
    Ok(s_line)
}

/// A MAF alignment record refer to https://genome.ucsc.edu/FAQ/FAQformat.html#format5
/// a pair of a-lines should be a align record
#[derive(Debug, PartialEq, Eq)]
pub struct MAFRecord {
    score: u64,
    slines: Vec<MAFSLine>,
}

// TODO: impl a derive macro for AlignRecord to cmp by target_start and target_name
impl PartialOrd<Self> for MAFRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.cmp(other).into()
    }
}

impl Ord for MAFRecord {
    fn cmp(&self, other: &Self) -> Ordering {
        let t1_name = self.target_name();
        let t1_start = self.target_start();
        let t2_name = other.target_name();
        let t2_start = other.target_start();
        if t1_name == t2_name {
            t1_start.cmp(&t2_start)
        } else {
            natord::compare(t1_name, t2_name)
        }
    }
}

/// impl Default trait for MAFRecord
impl Default for MAFRecord {
    fn default() -> Self {
        MAFRecord {
            score: 255,
            slines: Vec::new(),
        }
    }
}

/// A MAF record iterator
/// two s-lines should be a record
pub struct MAFRecords<'a, R: Read+Send> {
    inner: &'a mut BufReader<R>,
}

/// impl Iterator trait for MAFRecords
impl<R: Read+Send> Iterator for MAFRecords<'_, R> {
    type Item = Result<MAFRecord, ParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        let score = 255;
        match self.inner.lines().next() {
            Some(Ok(line)) => {
                if !line.starts_with('s') {
                    self.next() // skip empty line
                } else {
                    // start read multi s-lines
                    let mut mafrecord = MAFRecord {
                        // init a maf-record
                        score,
                        slines: Vec::new(),
                    };
                    let sline = match sline_from_string(line) {
                        Ok(sline) => sline,
                        Err(e) => return Some(Err(e)),
                    };
                    mafrecord.slines.push(sline); // push first s-line
                                                  // start read next sequential s-lines
                    for line in self.inner.lines() {
                        match line {
                            Ok(line) => {
                                if line.starts_with('s') {
                                    let sline = match sline_from_string(line) {
                                        Ok(sline) => sline,
                                        Err(e) => return Some(Err(e)),
                                    };
                                    mafrecord.slines.push(sline);
                                } else {
                                    // if s-line is over, break
                                    break;
                                }
                            }
                            _ => {
                                // if line is empty, break
                                break;
                            }
                        }
                    }
                    Some(Ok(mafrecord))
                }
            }
            _ => None, // if line is empty, iterator over
        }
    }
}

/// impl AlignRecord Trait for PafRecord
impl AlignRecord for MAFRecord {
    fn query_name(&self) -> &str {
        self.slines[1].name.as_str()
    }

    fn query_length(&self) -> u64 {
        self.slines[1].size
    }

    fn query_start(&self) -> u64 {
        self.slines[1].start
    }

    fn query_end(&self) -> u64 {
        self.slines[1].start + self.slines[1].align_size
    }

    fn query_strand(&self) -> Strand {
        self.slines[1].strand
    }

    fn target_name(&self) -> &str {
        self.slines[0].name.as_str()
    }

    fn target_length(&self) -> u64 {
        self.slines[0].size
    }

    fn target_start(&self) -> u64 {
        self.slines[0].start
    }

    fn target_end(&self) -> u64 {
        self.slines[0].start + self.slines[0].align_size
    }

    fn target_strand(&self) -> Strand {
        self.slines[0].strand
    }

    fn get_cigar_string(&self) -> String {
        parse_maf_seq_to_cigar(self, false).cigar_string
    }

    fn convert2paf(&self) -> PafRecord {
        let cigar = parse_maf_seq_to_cigar(self, false);
        let cigar_string = String::from("cg:Z:") + &cigar.cigar_string;
        let matches = cigar.match_count as u64;
        let block_length =
            (cigar.match_count + cigar.mismatch_count + cigar.ins_count + cigar.del_count) as u64;
        let edit_dist = cigar.mismatch_count + cigar.ins_count + cigar.del_count;
        let nm_tag = String::from("NM:i:") + &*edit_dist.to_string();

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
            tags: vec![nm_tag, cigar_string],
        }
    }

    fn convert2bam(&self, name_id_map: &HashMap<&str, u64>) -> SamRecord {
        // init a bam record
        let mut bamrec = SamRecord::default();

        // set bam record query name
        let q_name = self.query_name();
        let q_name: ReadName = q_name.parse().unwrap(); // TODO: handle parse error
        *bamrec.read_name_mut() = Some(q_name);

        // set bam record flags: it always in empty in whole genome alignment
        *bamrec.flags_mut() = Flags::empty();

        // set bam record reference sequence id ref to name-id-map
        let t_name = self.target_name();
        let t_id = *name_id_map.get(t_name).unwrap();
        *bamrec.reference_sequence_id_mut() = Some(t_id as usize);

        // set bam record position
        let t_start = self.target_start() + 1; // 0-based to 1-based
        *bamrec.alignment_start_mut() = Position::new(t_start as usize);

        // set bam record cigar
        let gen_cigar = parse_maf_seq_to_cigar(self, false);
        let cigar: Cigar = gen_cigar.cigar_string.parse().unwrap();
        *bamrec.cigar_mut() = cigar;

        // set bam record sequence
        let mut q_seq_ref = self.query_seq().to_string();
        q_seq_ref.retain(|c| c != '-'); // remove gap;should be UPPER?
        let q_seq: Sequence = q_seq_ref.parse().unwrap();
        *bamrec.sequence_mut() = q_seq;

        // set bam record tags
        let edit_dist = gen_cigar.mismatch_count + gen_cigar.ins_count + gen_cigar.del_count;
        let nm_tag = String::from("NM:i:") + &*edit_dist.to_string();
        let tag: Data = nm_tag.parse().unwrap();
        *bamrec.data_mut() = tag;

        // return bam record
        bamrec
    }

    fn query_seq(&self) -> &str {
        &self.slines[1].seq
    }

    fn target_seq(&self) -> &str {
        &self.slines[0].seq
    }
}
