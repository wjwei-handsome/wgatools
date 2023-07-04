//! The error kinds when process whole genome alignments(wga)

use std::fmt::Formatter;
use std::{fmt, io};

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
}

/// Represents what a file and it's format when error occurs
#[derive(Debug, PartialEq, Clone)]
pub struct FileInfo {
    /// The file name
    pub name: String,
    /// The file format
    pub format: FileFormat,
}

/// The error kinds when parse file/stream
#[derive(Debug, PartialEq, Clone)]
pub enum ParseErrorKind {
    /// Error when file or stream input/output
    Io,
    /// Empty file/stream
    Empty,
    /// Serde error when deserialize
    Serde,
    /// CIGAR parse error
    ParseCigar,
}

///The parse error returns
#[derive(Debug, PartialEq, Clone)]
pub struct ParseError {
    /// A description of what went wrong
    pub msg: String,
    /// The type of error that occurred
    pub kind: ParseErrorKind,
    /// The file and it's format when error occurs
    pub file_info: FileInfo,
}

impl ParseError {
    pub fn new_parse_cigar_err(format: FileFormat) -> Self {
        ParseError {
            msg: String::from("Parse cigar error, please check the format."),
            kind: ParseErrorKind::ParseCigar,
            file_info: FileInfo {
                name: String::from(""),
                format,
            },
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self.kind {
            ParseErrorKind::Io => {
                write!(
                    f,
                    "IO error: {}, check path: {}",
                    self.msg, self.file_info.name
                )
            },
            ParseErrorKind::Empty => {
                write!(
                    f,
                    "Empty input: {}, check path: {}",
                    self.msg, self.file_info.name
                )
            },
            _ => {
                write!(
                    f,
                    "Parse error: {}, check path: {}",
                    self.msg, self.file_info.name
                )
            },
        }
    }
}

impl From<io::Error> for ParseError {
    fn from(err: io::Error) -> ParseError {
        ParseError {
            msg: err.to_string(),
            kind: ParseErrorKind::Io,
            file_info: FileInfo {
                name: String::new(),
                format: FileFormat::Unknown,
            },
        }
    }
}

impl From<csv::Error> for ParseError {
    fn from(err: csv::Error) -> ParseError {
        ParseError {
            msg: err.to_string(),
            kind: ParseErrorKind::Serde,
            file_info: FileInfo {
                name: String::new(),
                format: FileFormat::Unknown,
            },
        }
    }
}
