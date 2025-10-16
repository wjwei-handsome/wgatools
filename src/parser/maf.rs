use crate::errors::{ParseMafErrKind, WGAError};
use crate::parser::cigar::parse_maf_seq_to_cigar;
use crate::parser::common::{AlignRecord, RecStat, Strand};
use crate::parser::paf::PafRecord;
use crate::utils::parse_str2u64;
use anyhow::anyhow;
use log::warn;
use regex::Regex;
use std::cmp::Ordering;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader, Read};

/// Parser for MAF file format
pub struct MAFReader<R: Read> {
    pub inner: BufReader<R>,
    pub header: String,
}

impl<R> MAFReader<R>
where
    R: Read + Send,
{
    /// Create a new MAF parser
    pub fn new(reader: R) -> Result<Self, WGAError> {
        let mut buf_reader = BufReader::new(reader);
        let mut header = String::new();
        buf_reader.read_line(&mut header)?;
        if !header.starts_with('#') {
            warn!("MAF Header is not start with `#`")
        }
        Ok(MAFReader {
            inner: buf_reader,
            header,
        })
    }

    /// Iterate over the records in the MAF file
    pub fn records(&mut self) -> MAFRecords<R> {
        MAFRecords {
            inner: self.inner.by_ref(),
        }
    }
}

impl MAFReader<File> {
    /// Create a new PAF parser from a file path
    pub fn from_path<P: AsRef<std::path::Path>>(path: P) -> Result<MAFReader<File>, WGAError> {
        match File::open(path.as_ref()) {
            Ok(file) => Ok(MAFReader::new(file)?),
            Err(_) => Err(WGAError::FileNotExist(path.as_ref().to_path_buf())),
        }
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
pub struct MAFSLine {
    pub mode: char,
    pub name: String,
    pub start: u64,
    pub align_size: u64,
    pub strand: Strand,
    pub size: u64,
    pub seq: String,
}

// impl mut for MAFSLine
impl MAFSLine {
    //        A--TCCGA
    // region:1  23456
    // colidx:01234567
    // this function will return the col index of the region pos
    fn get_col_coord(&self, pos: u64) -> u64 {
        let mut region_pos = 0;

        for (i, c) in self.seq.chars().enumerate() {
            if c != '-' {
                // if find the pos, return the col index
                if region_pos == pos {
                    return i as u64;
                }
                region_pos += 1;
            }
        }
        // if not find the pos, return the col index of the last char
        self.seq.len() as u64
    }

    // fn get_col_coord(&self, pos: u64) -> u64 {
    //     let mut col_coord = 0;
    //     let mut flag = 0;
    //     // skip '-'
    //     for (i, c) in self.seq.chars().enumerate() {
    //         if c == '-' {
    //             continue;
    //         } else {
    //             flag += 1;
    //             if flag == pos + 1 {
    //                 col_coord = i as u64;
    //                 break;
    //             }
    //         }
    //     }
    //     if col_coord == 0 {
    //         // if col_coord is 0, it means pos is the last position
    //         col_coord = self.seq.len() as u64;
    //     }

    //     col_coord
    // }

    pub fn set_start(&mut self, start: u64) {
        self.start = start;
    }
    pub fn set_align_size(&mut self, align_size: u64) {
        self.align_size = align_size;
    }
    pub fn set_strand(&mut self, strand: Strand) {
        self.strand = strand;
    }
    pub fn set_size(&mut self, size: u64) {
        self.size = size;
    }
    pub fn set_name(&mut self, name: String) {
        self.name = name;
    }
}

// main parse function for s-line
fn parse_sline(line: String) -> Result<MAFSLine, WGAError> {
    let mut iter = line.split_whitespace();
    let mode = match iter.next() {
        Some(mode) => mode
            .chars()
            .next()
            .ok_or(WGAError::Other(anyhow!("mode is empty"))),
        None => {
            return Err(WGAError::ParseMaf(ParseMafErrKind::FiledMissing(
                "mode".to_string(),
            )))
        }
    }?;
    let name = match iter.next() {
        Some(name) => name.to_string(),
        None => {
            return Err(WGAError::ParseMaf(ParseMafErrKind::FiledMissing(
                "name".to_string(),
            )))
        }
    };
    let start = match iter.next() {
        Some(start) => parse_str2u64(start)?,
        None => {
            return Err(WGAError::ParseMaf(ParseMafErrKind::FiledMissing(
                "start".to_string(),
            )))
        }
    };
    let align_size = match iter.next() {
        Some(align_size) => parse_str2u64(align_size)?,
        None => {
            return Err(WGAError::ParseMaf(ParseMafErrKind::FiledMissing(
                "align_size".to_string(),
            )))
        }
    };
    let strand = match iter.next() {
        Some(strand) => strand.parse::<Strand>()?,
        None => {
            return Err(WGAError::ParseMaf(ParseMafErrKind::FiledMissing(
                "strand".to_string(),
            )))
        }
    };
    let size = match iter.next() {
        Some(size) => parse_str2u64(size)?,
        None => {
            return Err(WGAError::ParseMaf(ParseMafErrKind::FiledMissing(
                "size".to_string(),
            )))
        }
    };
    let seq = match iter.next() {
        Some(seq) => seq.to_string(),
        None => {
            return Err(WGAError::ParseMaf(ParseMafErrKind::FiledMissing(
                "seq".to_string(),
            )))
        }
    };
    if iter.next().is_some() {
        return Err(WGAError::ParseMaf(ParseMafErrKind::SurplusField));
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

/// A MAF alignment record refer to https://genome.ucsc.edu/FAQ/FAQformat.html#format5
/// a pair of a-lines should be a align record
#[derive(Debug, PartialEq, Eq)]
pub struct MAFRecord {
    pub score: u64,
    pub slines: Vec<MAFSLine>,
    pub query_idx: usize,
}

impl MAFRecord {
    pub fn slice_block(&mut self, cut_start: u64, cut_end: u64, ord: usize) {
        let sline = &mut self.slines[ord];

        let cut_start_index = cut_start - sline.start;
        let cut_end_index = cut_end - sline.start;

        sline.set_start(cut_start);
        sline.set_align_size(cut_end - cut_start);

        let start_coord = sline.get_col_coord(cut_start_index);
        let end_coord = sline.get_col_coord(cut_end_index);
        sline.seq = sline.seq[start_coord as usize..end_coord as usize].to_string();

        let mut sline_idx_vec = (0..self.slines.len()).collect::<Vec<usize>>();
        sline_idx_vec.remove(ord);
        for sline in sline_idx_vec.iter() {
            let sline = &mut self.slines[*sline];
            let new_s_start = sline.start + cut_start_index;
            sline.set_start(new_s_start);
            let new_seq = sline.seq[start_coord as usize..end_coord as usize].to_string();
            let pre_align_size = end_coord - start_coord;
            let gap_size = new_seq.matches('-').count() as u64;
            sline.set_align_size(pre_align_size - gap_size);
            sline.seq = new_seq;
        }
    }

    pub fn rename(&mut self, prefixs: &[&str]) -> Result<(), WGAError> {
        // check prefixs length and slines length
        if prefixs.len() != self.slines.len() {
            return Err(WGAError::SLineCountNotMatch);
        }
        for (order, sline) in self.slines.iter_mut().enumerate() {
            let prefix = prefixs[order];
            let new_name = format!("{}{}", prefix, sline.name);
            sline.set_name(new_name);
        }
        Ok(())
    }

    pub fn get_query_idx_by_name(&self, query_name: &str) -> Option<usize> {
        self.slines.iter().position(|x| x.name == query_name)
    }

    pub fn get_query_idx_by_regex(&self, query_regex: &Regex) -> Option<usize> {
        self.slines
            .iter()
            .position(|x| query_regex.is_match(&x.name))
    }

    pub fn set_query_idx(&mut self, query_idx: usize) {
        self.query_idx = query_idx;
    }

    pub fn set_query_idx_by_name(&mut self, query_name: &str) -> Result<(), WGAError> {
        match self.get_query_idx_by_name(query_name) {
            Some(idx) => {
                self.query_idx = idx;
                Ok(())
            }
            None => Err(WGAError::QueryNameNotFound(query_name.to_string())),
        }
    }

    pub fn set_query_idx_by_regex(&mut self, query_regex: &Regex) -> Result<(), WGAError> {
        match self.get_query_idx_by_regex(query_regex) {
            Some(idx) => {
                self.query_idx = idx;
                Ok(())
            }
            None => Err(WGAError::QueryNameNotFound(query_regex.to_string())),
        }
    }

    // pub fn query_name_idx(&self, idx: usize) -> &str {
    //     self.slines[idx].name.as_str()
    // }

    // pub fn query_length_idx(&self, idx: usize) -> u64 {
    //     self.slines[idx].size
    // }

    // pub fn query_strand_idx(&self, idx: usize) -> Strand {
    //     self.slines[idx].strand
    // }

    // pub fn query_start_idx(&self, idx: usize) -> u64 {
    //     match self.query_strand_idx(idx) {
    //         Strand::Positive => self.slines[idx].start,
    //         Strand::Negative => {
    //             self.slines[idx].size - self.slines[idx].start - self.slines[idx].align_size
    //         }
    //     }
    // }

    // pub fn query_end_idx(&self, idx: usize) -> u64 {
    //     match self.query_strand_idx(idx) {
    //         Strand::Positive => self.slines[idx].start + self.slines[idx].align_size,
    //         Strand::Negative => self.slines[idx].size - self.slines[idx].start,
    //     }
    // }

    // pub fn query_seq_idx(&self, idx: usize) -> &str {
    //     &self.slines[idx].seq
    // }
}

// impl PartialEq for MAFRecord
impl PartialOrd<Self> for MAFRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

// impl Ord for MAFRecord
impl Ord for MAFRecord {
    // natural order
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
            query_idx: 1,
        }
    }
}

/// A MAF record iterator
/// two s-lines should be a record
pub struct MAFRecords<'a, R: Read + Send> {
    inner: &'a mut BufReader<R>,
}

/// impl Iterator trait for MAFRecords
impl<R: Read + Send> Iterator for MAFRecords<'_, R> {
    type Item = Result<MAFRecord, WGAError>;

    fn next(&mut self) -> Option<Self::Item> {
        let score = 255;
        match self.inner.lines().next() {
            Some(Ok(line)) => {
                if !line.starts_with('s') {
                    self.next() // skip empty line
                } else {
                    // start read multi s-lines
                    // init a maf-record
                    let mut mafrecord = MAFRecord {
                        score,
                        slines: Vec::new(),
                        query_idx: 1,
                    };
                    let sline = match parse_sline(line) {
                        Ok(sline) => sline,
                        // if catch error, return error
                        Err(e) => return Some(Err(e)),
                    };
                    mafrecord.slines.push(sline); // push first s-line
                                                  // start read next sequential s-lines
                    for line in self.inner.lines() {
                        match line {
                            Ok(line) => {
                                if line.starts_with('s') {
                                    let sline = match parse_sline(line) {
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
        self.slines[self.query_idx].name.as_str()
    }

    fn query_length(&self) -> u64 {
        self.slines[self.query_idx].size
    }

    fn query_start(&self) -> u64 {
        // self.slines[1].start
        let i = self.query_idx;
        match self.query_strand() {
            Strand::Positive => self.slines[i].start,
            Strand::Negative => {
                self.slines[i].size - self.slines[i].start - self.slines[i].align_size
            }
        }
    }

    fn query_end(&self) -> u64 {
        let i = self.query_idx;
        match self.query_strand() {
            Strand::Positive => self.slines[i].start + self.slines[i].align_size,
            Strand::Negative => self.slines[i].size - self.slines[i].start,
        }
    }

    fn query_strand(&self) -> Strand {
        self.slines[self.query_idx].strand
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

    fn target_align_size(&self) -> u64 {
        self.slines[0].align_size
    }

    fn get_cigar_string(&self) -> Result<String, WGAError> {
        Ok(parse_maf_seq_to_cigar(self, false).cigar_string)
    }

    fn convert2paf(&mut self, query_name: Option<&str>) -> Result<PafRecord, WGAError> {
        match query_name {
            Some(qname) => {
                self.set_query_idx_by_name(qname)?;
            }
            None => {
                // do nothing
            }
        }
        let cigar = parse_maf_seq_to_cigar(self, false);
        let cigar_string = String::from("cg:Z:") + &cigar.cigar_string;
        let matches = cigar.match_count as u64;
        let block_length = (cigar.match_count
            + cigar.mismatch_count
            + cigar.ins_count
            + cigar.inv_ins_count
            + cigar.del_count
            + cigar.inv_del_count) as u64;
        let edit_dist = block_length - matches;
        let nm_tag = String::from("NM:i:") + &*edit_dist.to_string();

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
            tags: vec![nm_tag, cigar_string],
        })
    }

    fn query_seq(&self) -> &str {
        // &self.slines[self.query_idx].seq
        if self.query_idx < self.slines.len() {
            &self.slines[self.query_idx].seq
        } else {
            ""
        }
    }

    fn target_seq(&self) -> &str {
        &self.slines[0].seq
    }

    fn get_stat(&self) -> Result<RecStat, WGAError> {
        // just convert cigar to stat
        let cigar = parse_maf_seq_to_cigar(self, false);
        Ok(RecStat::from(cigar))
    }
}

/// A MAF Writer
pub struct MAFWriter<W>
where
    W: Write,
{
    inner: W,
}

impl<W> MAFWriter<W>
where
    W: Write,
{
    /// Create a new MAF writer
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// write header
    pub fn write_header(&mut self, header: String) -> Result<(), WGAError> {
        writeln!(self.inner, "{}", header)?;
        Ok(())
    }

    /// write records
    pub fn write_record(&mut self, record: &MAFRecord) -> Result<(), WGAError> {
        // write a-line
        let a_line = format!("a score={}\n", record.score);
        write!(self.inner, "{}", a_line)?;
        for sline in record.slines.iter() {
            // write s-line
            let s_line = format!(
                "s\t{}\t{}\t{}\t{}\t{}\t{}",
                sline.name, sline.start, sline.align_size, sline.strand, sline.size, sline.seq
            );
            writeln!(self.inner, "{}", s_line)?;
        }
        // write a empty line
        writeln!(self.inner)?;
        Ok(())
    }
}
