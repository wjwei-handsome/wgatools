//! The error kinds when process whole genome alignments(wga)

use crate::parser::common::FileFormat;
use std::error::Error;
use std::fmt::Formatter;
use std::num::ParseIntError;
use std::{fmt, io};

/// Represents what a file and it's format when error occurs
#[derive(Debug, PartialEq, Clone)]
pub struct FileInfo {
    /// The file name
    pub name: String,
    /// The file format
    pub format: Option<FileFormat>,
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
    /// MAF s-line parse error
    ParseSLine,
    /// Parse Int error
    ParseInt,
}

impl fmt::Display for ParseErrorKind {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            ParseErrorKind::Io => write!(f, "IO"),
            ParseErrorKind::Empty => write!(f, "Empty"),
            ParseErrorKind::Serde => write!(f, "Serde"),
            ParseErrorKind::ParseCigar => write!(f, "ParseCigar"),
            ParseErrorKind::ParseSLine => write!(f, "ParseSLine"),
            ParseErrorKind::ParseInt => write!(f, "ParseInt"),
        }
    }
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
                format: Some(format),
            },
        }
    }

    pub fn new_parse_int_err(input: &str) -> Self {
        ParseError {
            msg: format!("Parse int error, please check the input: {}", input),
            kind: ParseErrorKind::ParseInt,
            file_info: FileInfo {
                name: String::from(""),
                format: None,
            },
        }
    }

    pub fn new_parse_maf_err(field: &str) -> Self {
        ParseError {
            msg: format!(
                "Parse maf error, please check this wrong field `{}`.",
                field
            ),
            kind: ParseErrorKind::ParseSLine,
            file_info: FileInfo {
                name: String::from(""),
                format: Some(FileFormat::Maf),
            },
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "ErrorKind: {} ErrorMsg: {}", self.kind, self.msg)
    }
}

impl From<io::Error> for ParseError {
    fn from(err: io::Error) -> ParseError {
        ParseError {
            msg: err.to_string(),
            kind: ParseErrorKind::Io,
            file_info: FileInfo {
                name: String::new(),
                format: None,
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
                format: None,
            },
        }
    }
}

impl From<ParseIntError> for ParseError {
    fn from(err: ParseIntError) -> ParseError {
        ParseError {
            msg: err.to_string(),
            kind: ParseErrorKind::ParseInt,
            file_info: FileInfo {
                name: String::new(),
                format: None,
            },
        }
    }
}

impl Error for ParseError {
    fn description(&self) -> &str {
        &self.msg
    }
}
