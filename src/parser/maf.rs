use crate::parser::common::{FileFormat, Strand};
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, Read};

/// Parser for MAF file format
pub struct MAFReader<R: io::Read> {
    inner: BufReader<R>,
    pub header: String,
}

impl<R> MAFReader<R>
where
    R: io::Read,
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
    pub fn convert(&mut self, _outputpath: &str, format: FileFormat) {
        match format {
            FileFormat::Chain => {}
            FileFormat::Blocks => {}
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
#[derive(Debug, PartialEq)]
struct MAFSLine {
    mode: char,
    name: String,
    start: u64,
    align_size: u64,
    strand: Strand,
    size: u64,
    seq: String,
}

impl From<String> for MAFSLine {
    fn from(value: String) -> Self {
        // if use nom to parse more fast?
        let mut iter = value.split_whitespace();
        let mode = iter.next().unwrap().chars().next().unwrap();
        let name = iter.next().unwrap().to_string();
        let start = iter.next().unwrap().parse::<u64>().unwrap();
        let align_size = iter.next().unwrap().parse::<u64>().unwrap();
        let strand = iter.next().unwrap().chars().next().unwrap();
        let size = iter.next().unwrap().parse::<u64>().unwrap();
        let seq = iter.next().unwrap().to_string();
        MAFSLine {
            mode,
            name,
            start,
            align_size,
            strand: Strand::from(strand),
            size,
            seq,
        }
    }
}

/// A MAF alignment record refer to https://genome.ucsc.edu/FAQ/FAQformat.html#format5
/// a pair of a-lines should be a align record
#[derive(Debug, PartialEq)]
pub struct MAFRecord {
    score: u64,
    slines: Vec<MAFSLine>,
}

/// A MAF record iterator
/// two s-lines should be a record
pub struct MAFRecords<'a, R: io::Read> {
    inner: &'a mut BufReader<R>,
}

impl Iterator for MAFRecords<'_, File> {
    type Item = MAFRecord;

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
                    mafrecord.slines.push(MAFSLine::from(line)); // push first s-line
                                                                 // start read next sequential s-lines
                    for line in self.inner.lines() {
                        match line {
                            Ok(line) => {
                                if line.starts_with('s') {
                                    mafrecord.slines.push(MAFSLine::from(line));
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
                    Some(mafrecord)
                }
            }
            _ => None, // if line is empty, iterator over
        }
    }
}
