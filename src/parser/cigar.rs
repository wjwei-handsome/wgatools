use crate::errors::WGAError;
use crate::parser::chain::{ChainDataLine, ChainRecord};
use crate::parser::common::{AlignRecord, Block};
use crate::utils::parse_str2u64;
use csv::Writer;
use itertools::Itertools;
use nom::bytes::complete::{tag, take_till, take_while};
use nom::error::Error;
use nom::multi::fold_many1;
use nom::{AsChar, IResult};
use std::io::Write;
use std::str;

// define Cigar struct
pub struct Cigar {
    pub cigar_string: String,
    pub match_count: usize,
    pub mismatch_count: usize,
    pub ins_event: usize,
    pub ins_count: usize,
    pub del_event: usize,
    pub del_count: usize,
    pub inv_ins_event: usize,
    pub inv_ins_count: usize,
    pub inv_del_event: usize,
    pub inv_del_count: usize,
    pub inv_event: usize,
}

/// CigarUnit is a atom operation in cigar string
#[derive(Debug)]
struct CigarUnit {
    op: char, // M, I, D, N, S, H, P, =, X operations
    len: u64, // length of the operation
}

struct CigarStrTuple<'a>(&'a str, &'a str);
/// Convert CigarStrTuple to CigarUnit
/// NOTE: We have carried out seemingly cumbersome operations here,
/// in fact, in order to better transmit errors due to
/// nom::fold_many1, we need to return IResult for FnMut
fn cst2cu(cst: CigarStrTuple) -> Result<CigarUnit, WGAError> {
    // just consume the first char of cst.0, will faster than cst.0.as_bytes()[0]
    let mut chars = cst.0.chars();
    let op = match chars.next() {
        Some(c) => c,
        None => return Err(WGAError::CigarOpInvalid(cst.0.to_string())),
    };
    // read cst.0 next again to check if it is a valid op
    if chars.next().is_some() {
        return Err(WGAError::CigarOpInvalid(cst.0.to_string()));
    }
    let len = parse_str2u64(cst.1)?; // actually it will never occur error
    Ok(CigarUnit { op, len })
}

/// Parse a cigar unit until the input is empty
fn parse_cigar_str_tuple(input: &str) -> IResult<&str, CigarStrTuple> {
    // if input is empty, return error to break infinite loop
    if input.is_empty() {
        return Err(nom::Err::Error(Error::new(
            input,
            nom::error::ErrorKind::Eof,
        )));
    }

    // let input = input.as_bytes();
    // take digits
    let (input, len) = take_while(AsChar::is_dec_digit)(input)?;
    // take alphabetic
    let (input, op) = take_till(AsChar::is_dec_digit)(input)?;

    Ok((input, CigarStrTuple(op, len)))
}

/// a phantom
fn null() -> Result<(), WGAError> {
    Ok(())
}

fn cigar_unit_block(
    op: char,
    count: u64,
    wtr: &mut Writer<&mut dyn Write>,
    block: &mut Block,
) -> Result<(), WGAError> {
    match op {
        'M' => {
            // move query&target
            block.query_end += count;
            block.target_end += count;

            // write and serialize
            wtr.serialize(&block)?;

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
        _ => return Err(WGAError::CigarOpInvalid(op.to_string())), // TODO: handle `H` for SAM
    };
    Ok(())
}

/// Parse cigar string of a AlignRecord[PafRecord, SamRecord] which includes cg:Z: tag and
/// write into a blocks file
/// - For PafRecord: cigar should only contains 'M,I,D'
/// - For SamRecord: cigar's first `[0-9]+H` should represent the query start
pub fn parse_cigar_to_blocks<T: AlignRecord>(
    rec: &T,
    wtr: &mut Writer<&mut dyn Write>,
) -> Result<(), WGAError> {
    // get cigar bytes and tags
    let cigar = rec.get_cigar_str()?;
    let (cigar, _tag) = tag("cg:Z:")(cigar)?;

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
    let (_, res) = fold_many1(
        parse_cigar_str_tuple,
        null,
        |res: Result<(), WGAError>, cigarunit| {
            if res.is_ok() {
                let cigarunit = cst2cu(cigarunit)?;
                cigar_unit_block(cigarunit.op, cigarunit.len, wtr, &mut block)?
            }
            res
        },
    )(cigar)?;
    res
}

/// Parse cigar string of a AlignRecord[PafRecord, SamRecord] which includes cg:Z: tag and
/// write into a chain file.
/// - For PafRecord: cigar should only contains 'M,I,D'
/// - For SamRecord: cigar's first `[0-9]+H` should represent the query start
pub fn parse_cigar_to_chain<T: AlignRecord>(
    rec: &T,
    wtr: &mut Box<dyn Write>,
    // ) -> Result<(&'a str, Result<(), WGAError>), WGAError> {
) -> Result<(), WGAError> {
    // get cigar bytes and tag
    let cigar = rec.get_cigar_str()?;
    let (cigar, _tag) = tag("cg:Z:")(cigar)?;

    // init a ChainDataLine filled 0
    let mut dataline = ChainDataLine {
        size: 0,
        query_diff: 0,
        target_diff: 0,
    };

    // try regex for cigar
    // let re = regex::Regex::new(r"([0-9]+)([MID=X])")?;
    // for cap in re.captures_iter(cigar) {
    //     let len = cap[1].parse::<u64>()?;
    //     let op = cap[2].chars().next()?;
    //     cigar_unit_chain(op, len, wtr, &mut dataline)?;
    // }
    // Unfortunately, it turns out that this is three times slower than nom

    // fold cigar bytes into many CigarUnits[#CigarUnit]
    let (_, res) = fold_many1(
        parse_cigar_str_tuple,
        null,
        |res: Result<(), WGAError>, cigarunit| {
            // Continue execution only when no error occurs
            if res.is_ok() {
                let cigarunit = cst2cu(cigarunit)?;
                cigar_unit_chain(cigarunit.op, cigarunit.len, wtr, &mut dataline)?;
            }
            res
        },
    )(cigar)?;

    // After all cigar units successfully write done, the last dataline.size should be wrote
    if res.is_ok() {
        wtr.write_all(format!("\n{}", dataline.size).as_bytes())?;
    }
    res
}

/// cigar category method -- extension
pub fn cigar_cat_ext(c1: &char, c2: &char) -> char {
    if c1 == c2 {
        '='
    } else if c1 == &'-' {
        'I'
    } else if c2 == &'-' {
        'D'
    } else {
        'X'
    }
}

/// cigar category method
pub fn cigar_cat(c1: &char, c2: &char) -> char {
    if c1 == c2 {
        'M'
    } else if c1 == &'-' {
        'I'
    } else if c2 == &'-' {
        'D'
    } else {
        'M'
    }
}

/// parse MAF two seqs into Cigar
pub fn parse_maf_seq_to_cigar<T: AlignRecord>(rec: &T, with_h: bool) -> Cigar {
    let mut cigar_string = String::new();
    let seq1_iter = rec.target_seq().chars();
    let seq2_iter = rec.query_seq().chars();
    let mut match_count = 0;
    let mut mismatch_count = 0;
    let mut ins_event = 0;
    let mut ins_count = 0;
    let mut del_event = 0;
    let mut del_count = 0;
    let mut inv_ins_event = 0;
    let mut inv_ins_count = 0;
    let mut inv_del_event = 0;
    let mut inv_del_count = 0;
    let mut inv_event = 0;
    let group_by_iter = seq1_iter
        .zip(seq2_iter)
        .group_by(|(c1, c2)| cigar_cat_ext(c1, c2));

    let begin = rec.query_start() as usize;
    let end = rec.query_length() - rec.query_end();
    if with_h {
        cigar_string.push_str(&begin.to_string());
        cigar_string.push('H');
    }

    let inv = match rec.query_strand() {
        crate::parser::common::Strand::Positive => false,
        crate::parser::common::Strand::Negative => {
            inv_event = 1;
            true
        }
    };

    for (k, g) in group_by_iter.into_iter() {
        let len = g.count();
        // 10=5X1D2I ==> 15M1D2I
        // (10,=),(5,X),(1D),(2I)
        match k {
            '=' => {
                match_count += len;
            }
            'I' => {
                if inv {
                    inv_ins_event += 1;
                    inv_ins_count += len;
                } else {
                    ins_event += 1;
                    ins_count += len;
                }
            }
            'D' => {
                if inv {
                    inv_del_event += 1;
                    inv_del_count += len;
                } else {
                    del_event += 1;
                    del_count += len;
                }
            }
            'X' => {
                mismatch_count += len;
            }
            _ => {}
        };
        cigar_string.push_str(&len.to_string());
        cigar_string.push(k);
    }

    if with_h {
        cigar_string.push_str(&end.to_string());
        cigar_string.push('H');
    }

    Cigar {
        cigar_string,
        match_count,
        mismatch_count,
        ins_event,
        ins_count,
        del_event,
        del_count,
        inv_ins_event,
        inv_ins_count,
        inv_del_event,
        inv_del_count,
        inv_event,
    }
}

/// parse MAF two seqs adn write into a chain file
pub fn parse_maf_seq_to_chain<T: AlignRecord>(
    rec: &T,
    wtr: &mut Box<dyn Write>,
) -> Result<(), WGAError> {
    let seq1_iter = rec.target_seq().chars();
    let seq2_iter = rec.query_seq().chars();
    let group_by_iter = seq1_iter
        .zip(seq2_iter)
        .group_by(|(c1, c2)| cigar_cat(c1, c2));

    // init a ChainDataLine filled 0
    let mut dataline = ChainDataLine {
        size: 0,
        query_diff: 0,
        target_diff: 0,
    };
    for (k, g) in group_by_iter.into_iter() {
        let len = g.count() as u64;
        cigar_unit_chain(k, len, wtr, &mut dataline)?;
    }
    // After all cigar units write done, the last dataline.size should be wrote
    wtr.write_all(format!("\n{}", dataline.size).as_bytes())?;
    Ok(())
}

fn cigar_unit_chain(
    op: char,
    count: u64,
    wtr: &mut Box<dyn Write>,
    dataline: &mut ChainDataLine,
) -> Result<(), WGAError> {
    match op {
        'M' | 'X' | '=' => {
            // will not write unless: [1. size == 0; 2. both no query&target diff]
            if (dataline.size != 0) && (dataline.target_diff + dataline.query_diff != 0) {
                wtr.write_all(format!("{}", dataline).as_bytes())?;
                dataline.size = 0;
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
        _ => return Err(WGAError::CigarOpInvalid(op.to_string())), // TODO: handle 'H' for SAM
    };
    Ok(())
}

fn cigar_unit_insert_seq(
    op: char,
    count: u64,
    current_offset: &mut u64,
    t_seq: &mut String,
    q_seq: &mut String,
) -> Result<(), WGAError> {
    match op {
        'M' | '=' | 'X' => {
            // do nothing but move offset
            *current_offset += count;
        }
        'I' => {
            // insert '-' into target seq
            let ins_str = "-".repeat(count as usize);
            t_seq.insert_str(*current_offset as usize, &ins_str);
            *current_offset += count;
        }
        'D' => {
            // insert '-' into query seq
            let del_str = "-".repeat(count as usize);
            q_seq.insert_str(*current_offset as usize, &del_str);
            *current_offset += count;
        }
        _ => return Err(WGAError::CigarOpInvalid(op.to_string())), // TODO: handle 'H' for SAM
    };
    Ok(())
}

/// Parse cigar to insert `-` in MAF sequences
pub fn parse_cigar_to_insert<T: AlignRecord>(
    rec: &T,
    t_seq: &mut String,
    q_seq: &mut String,
) -> Result<(), WGAError> {
    // get cigar bytes and tag
    let cigar = rec.get_cigar_str()?;
    let (cigar, _tag) = tag("cg:Z:")(cigar)?;

    // fold cigar bytes into many CigarUnits[#CigarUnit]
    let mut current_offset = 0;
    let (_, res) = fold_many1(
        parse_cigar_str_tuple,
        null,
        |res: Result<(), WGAError>, cigarunit| {
            if res.is_ok() {
                let cigarunit = cst2cu(cigarunit)?;
                cigar_unit_insert_seq(
                    cigarunit.op,
                    cigarunit.len,
                    &mut current_offset,
                    t_seq,
                    q_seq,
                )?
            }
            res
        },
    )(cigar)?;
    res
}

/// parse ChainRecord into Cigar
pub fn parse_chain_to_cigar(rec: &ChainRecord, _with_h: bool) -> Cigar {
    let mut cigar_string = String::new();
    let mut match_count = 0;
    let mismatch_count = 0;
    let mut ins_event = 0;
    let mut ins_count = 0;
    let mut del_event = 0;
    let mut del_count = 0;
    let mut inv_ins_event = 0;
    let mut inv_ins_count = 0;
    let mut inv_del_event = 0;
    let mut inv_del_count = 0;
    let mut inv_event = 0;

    let inv = match rec.query_strand() {
        crate::parser::common::Strand::Positive => false,
        crate::parser::common::Strand::Negative => {
            inv_event = 1;
            true
        }
    };

    for dataline in &rec.lines {
        let match_len = dataline.size;
        let ins_len = dataline.target_diff;
        let del_len = dataline.query_diff;
        cigar_string.push_str(&match_len.to_string());
        cigar_string.push('M');
        match_count += match_len as usize;
        match ins_len {
            0 => {}
            _ => {
                cigar_string.push_str(&ins_len.to_string());
                cigar_string.push('I');
                if inv {
                    inv_ins_event += 1;
                    inv_ins_count += ins_len as usize;
                } else {
                    ins_event += 1;
                    ins_count += ins_len as usize;
                }
            }
        };
        match del_len {
            0 => {}
            _ => {
                cigar_string.push_str(&del_len.to_string());
                cigar_string.push('D');
                if inv {
                    inv_del_event += 1;
                    inv_del_count += del_len as usize;
                } else {
                    del_event += 1;
                    del_count += del_len as usize;
                }
            }
        }
    }
    Cigar {
        cigar_string,
        match_count,
        mismatch_count,
        ins_event,
        ins_count,
        del_event,
        del_count,
        inv_ins_event,
        inv_ins_count,
        inv_del_event,
        inv_del_count,
        inv_event,
    }
}

/// Parse CIGAR to Cigar struct
pub fn parse_paf_to_cigar<T: AlignRecord>(rec: &T) -> Result<Cigar, WGAError> {
    let cigar_string = String::new();
    let mut match_count = 0;
    let mut mismatch_count = 0;
    let mut ins_event = 0;
    let mut ins_count = 0;
    let mut del_event = 0;
    let mut del_count = 0;
    let mut inv_ins_event = 0;
    let mut inv_ins_count = 0;
    let mut inv_del_event = 0;
    let mut inv_del_count = 0;
    let mut inv_event = 0;

    let inv = match rec.query_strand() {
        crate::parser::common::Strand::Positive => false,
        crate::parser::common::Strand::Negative => {
            inv_event = 1;
            true
        }
    };

    let cigar = rec.get_cigar_str()?;
    let (cigar, _tag) = tag("cg:Z:")(cigar)?;
    let mut current_offset = 0;
    let (_, res) = fold_many1(
        parse_cigar_str_tuple,
        null,
        |res: Result<(), WGAError>, cigarunit| {
            if res.is_ok() {
                let cigarunit = cst2cu(cigarunit)?;
                match cigarunit.op {
                    'M' | '=' => {
                        match_count += cigarunit.len as usize;
                    }
                    'X' => {
                        mismatch_count += cigarunit.len as usize;
                    }
                    'I' => {
                        if inv {
                            inv_ins_event += 1;
                            inv_ins_count += cigarunit.len as usize;
                        } else {
                            ins_event += 1;
                            ins_count += cigarunit.len as usize;
                        }
                    }
                    'D' => {
                        if inv {
                            inv_del_event += 1;
                            inv_del_count += cigarunit.len as usize;
                        } else {
                            del_event += 1;
                            del_count += cigarunit.len as usize;
                        }
                    }
                    _ => return Err(WGAError::CigarOpInvalid(cigarunit.op.to_string())), // TODO: handle 'H' for SAM
                };
                current_offset += cigarunit.len;
            }
            res
        },
    )(cigar)?;
    res?;
    Ok(Cigar {
        cigar_string,
        match_count,
        mismatch_count,
        ins_event,
        ins_count,
        del_event,
        del_count,
        inv_ins_event,
        inv_ins_count,
        inv_del_event,
        inv_del_count,
        inv_event,
    })
}

/// Parse CIGAR to Cigar struct and stat cov
pub fn update_cov_vec(cov_vec: &mut Vec<usize>, cigar: &str, start: usize) -> Result<(), WGAError> {
    let (cigar, _tag) = tag("cg:Z:")(cigar)?;
    let mut pos = start;
    let (_, res) = fold_many1(
        parse_cigar_str_tuple,
        null,
        |res: Result<(), WGAError>, cigarunit| {
            if res.is_ok() {
                let cigarunit = cst2cu(cigarunit)?;
                let length = cigarunit.len as usize;
                match cigarunit.op {
                    'M' | '=' => {
                        for i in pos..(pos + length) {
                            if i < cov_vec.len() {
                                cov_vec[i] += 1;
                            }
                        }
                        pos += length;
                    }
                    'I' | 'S' => {}

                    _ => {
                        pos += length;
                    }
                };
            }
            res
        },
    )(cigar)?;
    res?;
    Ok(())
}
