//! The error kinds when process whole genome alignments(wga)

use crate::tools::mafextra::GenomeRegion;
use thiserror::Error;

// define Error types
#[derive(Error, Debug)]
pub enum WGAError {
    // IO error
    #[error("IO error:{0}")]
    Io(#[from] std::io::Error),
    // File path not exist
    #[error("File path `{0}` not exist")]
    FileNotExist(std::path::PathBuf),
    // Parse MAF Error
    #[error("Parse MAF error by: {0}")]
    ParseMaf(ParseMafErrKind),
    #[error("CSV deserialize error by: {0}")]
    CsvDeserialize(#[from] csv::Error),
    #[error("Empty stdin, please add `-h` for help")]
    EmptyStdin,
    #[error("File `{0}` already exists, please add `-r` to rewrite it.")]
    FileReWrite(String),
    #[error("json dese error")]
    SerdeDeserialize(#[from] serde_json::Error),
    #[error("ThreadPoolBuildError error")]
    ThreadPoolBuildError(#[from] rayon::ThreadPoolBuildError),
    #[error("Empty record")]
    EmptyRecord,
    #[error("regions or region_file must be specified")]
    EmptyRegion,
    #[error("Stdin not allowed here")]
    StdinNotAllowed,
    #[error("Parse Genome Region Error By: {0}")]
    ParseGenomeRegion(ParseGenomeRegionErrKind),
    #[error("Failed region: {0}")]
    FailedRegion(GenomeRegion),
    #[error("Duplicate name `{0}` in a record not allowed, please check or use `rename`")]
    DuplicateName(String),
    #[error("Format {0} Parse Error by rust::nom, please check")]
    NomErr(#[from] nom::error::Error<String>),
    #[error("Parse Chain Error By: {0}")]
    ParseChain(ParseChainErrKind),
    #[error("Parse Strand `{0}` Error")]
    ParseStrand(String),
    #[error("Parse `{0}` Into Integer Error")]
    #[from = "ParseIntError"]
    ParseIntError(String),
    #[error("Parse `{0}` Into Float Error")]
    #[from = "ParseFloatError"]
    ParseFloatError(String),
    #[error("CIGAR start tag not found")]
    CigarTagNotFound,
    #[error("CIGAR OP `{0}` invalid")]
    CigarOpInvalid(String),
    #[error("noodles-sam parse error {0}")]
    NoodlesSamParseError(#[from] noodles::sam::record::reference_sequence_name::ParseError),
    #[error("noodlesp-sam try into num parse error {0}")]
    TryIntoNum(#[from] std::num::TryFromIntError),
    #[error("noodlesp-sam read name parse error {0}")]
    ReadNameParseError(#[from] noodles::sam::record::read_name::ParseError),
    #[error("HTS library error by {0}")]
    HtsLibError(#[from] rust_htslib::errors::Error),
    #[error("Unexcepted Regex Error by: {0}")]
    UnexceptedRegexError(String),
    #[error("Regex build Error")]
    RegexBuildError(#[from] regex::Error),
    #[error("Invalid Base: `{0}`")]
    InvalidBase(String),
    #[error("Ah-oh! NOT IMPLEMENTED :(")]
    NotImplemented,
    #[error("S-line count not match")]
    SLineCountNotMatch,
    // Other error
    #[error(transparent)]
    Other(#[from] anyhow::Error),
}

impl From<nom::Err<nom::error::Error<&str>>> for WGAError {
    fn from(value: nom::Err<nom::error::Error<&str>>) -> Self {
        match value {
            nom::Err::Error(e) => {
                WGAError::NomErr(nom::error::Error::new(e.input[..10].to_string(), e.code))
            }
            _ => WGAError::Other(anyhow::anyhow!("Other nom Error")),
        }
    }
}

#[derive(Error, Debug)]

pub enum ParseMafErrKind {
    #[error("S-line Filed `{0}` Missing")]
    FiledMissing(String),
    #[error("Surplus Filed > 7")]
    SurplusField,
}

#[derive(Error, Debug)]
pub enum ParseChainErrKind {
    #[error("Chain Line Field `{0}` Missing")]
    FiledMissing(String),
}

#[derive(Error, Debug)]
pub enum ParseGenomeRegionErrKind {
    #[error("Region `{0}` is match the format of `chr:start-end`")]
    FormatNotMatch(String),
    #[error("Start `{0}` is larger than end `{1}`")]
    StartGTEnd(u64, u64),
}
