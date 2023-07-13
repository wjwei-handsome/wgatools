use crate::parser::chain::ChainDataLine;
use crate::parser::common::{AlignRecord, Block};
use csv::Writer;
use nom::bytes::complete::{tag, take_while};
use nom::character::{is_alphabetic, is_digit};
use nom::error::Error;
use nom::multi::fold_many0;
use nom::IResult;
use std::io::Write;
use std::str;
use itertools::Itertools;

/// CigarUnit is a atom operation in cigar string
#[derive(Debug)]
pub struct CigarUnit {
    op: char, // M, I, D, N, S, H, P, =, X operations
    len: u64, // length of the operation
}

/// Parse a cigar unit until the input is empty
fn parse_cigar_unit(input: &[u8]) -> IResult<&[u8], CigarUnit> {
    // if input is empty, return error to break infinite loop
    if input.is_empty() {
        return Err(nom::Err::Error(Error::new(
            input,
            nom::error::ErrorKind::Eof,
        )));
    }

    // take digits
    let (input, len) = take_while(is_digit)(input)?;
    // take alphabetic
    let (input, op) = take_while(is_alphabetic)(input)?;
    // convert len to u64
    // TODO: handle parse error
    let len = str::from_utf8(len).unwrap().parse::<u64>().unwrap();

    let cigar_unit = CigarUnit {
        // TODO: op.len() always === 1, if need a parse err? or pass to match op {}?
        op: char::from(op[0]),
        len,
    };
    Ok((input, cigar_unit))
}

/// a phantom
fn null() {}

/// Parse cigar string of a AlignRecord[PafRecord, SamRecord] which includes cg:Z: tag and
/// write into a blocks file
/// - For PafRecord: cigar should only contains 'M,I,D'
/// - For SamRecord: cigar's first `[0-9]+H` should represent the query start
pub fn parse_cigar_to_blocks<'a, T: AlignRecord>(
    rec: &'a T,
    wtr: &mut Writer<Box<dyn Write>>,
) -> IResult<&'a [u8], ()> {
    // get cigar bytes and tags
    let cigar = rec.get_cigar_bytes();
    let (cigar, _tag) = tag(b"cg:Z:")(cigar)?;

    // init a original block
    let mut block = Block {
        query_name: rec.query_name(),
        query_start: rec.query_start(),
        query_end: rec.query_start(),
        target_name: rec.target_name(),
        target_start: rec.target_start(),
        target_end: rec.target_start(),
        strand: rec.query_strand(),
    };

    // fold cigar bytes into many CigarUnits[#CigarUnit]
    let (rest, res) = fold_many0(parse_cigar_unit, null, |(), cigarunit| {
        // get operate counts
        let count = cigarunit.len;
        // match operations
        match cigarunit.op {
            'M' => {
                // move query&target
                block.query_end += count;
                block.target_end += count;

                // write and serialize
                wtr.serialize(block).unwrap(); // TODO: handle IO error

                // sync query&target start with end
                block.query_start = block.query_end;
                block.target_start = block.target_end;
            }
            'I' => {
                // only move query start&end
                block.query_end += count;
                block.query_start += count;
            }
            'D' => {
                // only move target start&end
                block.target_end += count;
                block.target_start += count;
            }
            _ => {} // TODO: handle `H` for SAM
        };
    })(cigar)?;

    Ok((rest, res))
}

/// Parse cigar string of a AlignRecord[PafRecord, SamRecord] which includes cg:Z: tag and
/// write into a chain file.
/// - For PafRecord: cigar should only contains 'M,I,D'
/// - For SamRecord: cigar's first `[0-9]+H` should represent the query start
pub fn parse_cigar_to_chain<'a, T: AlignRecord>(
    rec: &'a T,
    wtr: &mut Box<dyn Write>,
) -> IResult<&'a [u8], ()> {
    // get cigar bytes and tag
    let cigar = rec.get_cigar_bytes();
    let (cigar, _tag) = tag(b"cg:Z:")(cigar)?;

    // init a ChainDataLine filled 0
    let mut dataline = ChainDataLine {
        size: 0,
        query_diff: 0,
        target_diff: 0,
    };

    // fold cigar bytes into many CigarUnits[#CigarUnit]
    let (rest, res) = fold_many0(parse_cigar_unit, null, |(), cigarunit| {
        // get operate counts
        let count = cigarunit.len;
        // match operations
        match cigarunit.op {
            'M' => {
                // will not write unless: [1. size == 0; 2. both no query&target diff]
                if (dataline.size != 0) && (dataline.target_diff + dataline.query_diff != 0) {
                    // TODO: handle IO error
                    wtr.write_all(format!("{}", dataline).as_bytes()).unwrap();
                };
                // accumulate size
                dataline.size += count;
                // init query&target diff
                dataline.target_diff = 0;
                dataline.query_diff = 0;
            }
            'I' => {
                // accumulate target diff for 'I'
                dataline.target_diff += count;
            }
            'D' => {
                // accumulate query diff for 'D'
                dataline.query_diff += count;
            }
            _ => {} // TODO: handle 'H' for SAM
        };
    })(cigar)?;

    // After all cigar units write done, the last dataline.size should be wrote
    wtr.write_all(format!("\n{}", dataline.size).as_bytes())
        .unwrap(); // TODO: handle IO error
    Ok((rest, res))
}

/// cigar category method
fn cigar_cat(c1: &char, c2: &char) -> &'static str {
    if c1 == c2 {
        "M"
    } else if c1 == &'-' {
        "I"
    } else if c2 == &'-' {
        "D"
    } else {
        "M"
    }
}

/// parse MAF two seqs into cigar string
pub fn parse_maf_seq_to_cigar<T: AlignRecord>(rec: &T) -> String {
    let mut cigar = String::new();
    let seq1_iter = rec.target_seq().chars();
    let seq2_iter = rec.query_seq().chars();
    seq1_iter
        .zip(seq2_iter)
        .group_by(|(c1, c2)| cigar_cat(c1, c2))
        .into_iter()
        .for_each(|(k, g)| {
            let len = g.count();
            cigar.push_str(&len.to_string());
            cigar.push_str(k);
        });
    cigar
}