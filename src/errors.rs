//! The error kinds when process whole genome alignments(wga)

use crate::tools::mafextra::GenomeRegion;
use thiserror::Error;

// define Error types
#[derive(Error, Debug)]
pub enum WGAError {
    // IO error
    #[error("IO Error:{0}")]
    Io(#[from] std::io::Error),
    // Parse MAF Error
    #[error("Parse MAF Error By: {0}")]
    ParseMaf(ParseMafErrKind),
    #[error("csv dese error")]
    CsvDeserialize(#[from] csv::Error),
    #[error("Empty stdin, please add `-h` for help")]
    EmptyStdin,
    #[error("File `{0}` already exists, please add `-r` to rewrite it.")]
    FileReWrite(String),
    #[error(transparent)]
    Other(#[from] anyhow::Error),
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
    #[error("failed region: {0}")]
    FailedRegion(GenomeRegion),
    #[error("Duplicate name `{0}` in a record not allowed, please check or use `modname`")]
    DuplicateName(String),
    #[error("Format {0} Parse Error, please check")]
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
    // #[error("Parse `{0}` Into Integer Error")]
    // #[from = "ParseIntError"]
    // ParseIntError(String),
}

#[derive(Error, Debug)]
pub enum ParseChainErrKind {
    #[error("Chain Line Field `{0}` Missing")]
    FiledMissing(String),
    // #[error("Parse `{0}` Into Integer Error")]
    // #[from = "ParseIntError"]
    // ParseIntError(String),
}

#[derive(Error, Debug)]
pub enum ParseGenomeRegionErrKind {
    #[error("Region `{0}` is match the format of `chr:start-end`")]
    FormatNotMatch(String),
    #[error("Start `{0}` is larger than end `{1}`")]
    StartGTEnd(u64, u64),
}
