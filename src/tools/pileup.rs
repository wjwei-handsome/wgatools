use crate::{
    errors::WGAError,
    parser::{
        common::AlignRecord,
        maf::{MAFReader, MAFRecord},
    },
};
use itertools::Itertools;
// use noodles::vcf;
// use noodles::vcf::{
//     header::{
//         record::value::{
//             map::{format::Type as fmttype, info::Type as infotype, Contig, Format, Info},
//             Map,
//         },
//         Number,
//     },
//     record::{
//         genotypes::keys::key as gtkey, info::field::key as infokey, Info as recinfo, Position,
//     },
//     Header, Record,
// };
// use rayon::prelude::*;
use std::{
    collections::{HashMap, HashSet},
    io::{Read, Write},
};

#[derive(Debug)]
struct Pileup {
    ref_name: String,
    ref_pos: u64,
    ref_base: char,
    alt: Alt,
    uid: String,
}

impl std::fmt::Display for Pileup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.alt.vartype {
            VarType::Ins | VarType::Del => write!(
                f,
                "{}\t{}\t{}\t{}{}\t{}",
                self.ref_name, self.ref_pos, self.ref_base, self.ref_base, self.alt, self.uid,
            ),
            VarType::Snp | VarType::Null => write!(
                f,
                "{}\t{}\t{}\t{}\t{}",
                self.ref_name, self.ref_pos, self.ref_base, self.alt, self.uid,
            ),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct Alt {
    len: u64,
    vartype: VarType,
    alt_base: String,
}

impl std::fmt::Display for Alt {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.vartype {
            VarType::Ins => write!(f, "+{}{}", self.len, self.alt_base),
            VarType::Del => write!(f, "-{}{}", self.len, self.alt_base),
            VarType::Snp | VarType::Null => write!(f, "{}", self.alt_base),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
enum VarType {
    Ins,
    Del,
    Snp,
    Null,
}

pub fn pileup<R: Read + Send>(
    mafreader: &mut MAFReader<R>,
    _writer: &mut dyn Write,
    all: bool,
) -> Result<(), WGAError> {
    let mut all_pileup_vecs = Vec::new();
    for rec in mafreader.records() {
        let rec = rec?;
        let pileup_vec = generate_pileup(&rec, all)?;
        all_pileup_vecs.push(pileup_vec);
    }
    let res = merge_pileup_vec(all_pileup_vecs)?;
    for i in res {
        println!("{:?}", i)
    }

    Ok(())
}

#[derive(Debug)]
struct MergedPileup {
    chro: String,
    pos: u64,
    ref_base: char,
    alts: Vec<Alt>,
    gts: HashMap<String, String>,
}

/// merge pileup vecs into hashmap
fn merge_pileup_vec(pileup_lists: Vec<Vec<Pileup>>) -> Result<Vec<MergedPileup>, WGAError> {
    let mut merged_map: HashMap<(String, u64), (char, HashSet<Alt>, HashMap<String, String>)> =
        HashMap::new();

    for list in pileup_lists {
        for pileup in list {
            let key = (pileup.ref_name.clone(), pileup.ref_pos);
            let entry =
                merged_map
                    .entry(key)
                    .or_insert((pileup.ref_base, HashSet::new(), HashMap::new()));

            entry.1.insert(pileup.alt.clone());

            let alt_index = entry.1.iter().position(|alt| *alt == pileup.alt).unwrap();
            let gt_value = format!("{}/{}", alt_index + 1, alt_index + 1);

            entry.2.insert(pileup.uid, gt_value);
        }
    }

    let mut merged_variants: Vec<MergedPileup> = merged_map
        .into_iter()
        .map(|((chro, pos), (ref_base, alts_set, gts))| MergedPileup {
            chro,
            pos,
            ref_base,
            alts: alts_set.into_iter().collect::<Vec<Alt>>(),
            gts,
        })
        .collect();

    merged_variants.sort_by_key(|mv| (mv.chro.clone(), mv.pos));
    Ok(merged_variants)
}

/// cigar category method -- extension
fn cigar_cat_ext(c1: &char, c2: &char) -> char {
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

/// generate pileup for a MAF record
fn generate_pileup(rec: &MAFRecord, all: bool) -> Result<Vec<Pileup>, WGAError> {
    let ref_name = rec.target_name();
    let ref_start = rec.target_start();
    let ref_end = rec.target_end();
    let mut ref_pos = ref_start;
    let mut ref_offset = 0;
    let query_name = rec.query_name();
    let query_start = rec.query_start();
    let query_end = rec.query_end();
    let uid = format!(
        "{}#{}#{}@{}#{}#{}",
        ref_name, ref_start, ref_end, query_name, query_start, query_end
    );

    let mut pileup_vec = Vec::new();
    let seq1_iter = rec.target_seq().chars();
    let seq2_iter = rec.query_seq().chars();
    let group_by_iter = seq1_iter
        .zip(seq2_iter)
        .group_by(|(c1, c2)| cigar_cat_ext(c1, c2));

    for (k, g) in group_by_iter.into_iter() {
        let g: Vec<(char, char)> = g.collect();
        let len = g.len() as u64;

        match k {
            '=' => {
                for _ in 0..len {
                    ref_pos += 1;
                    ref_offset += 1;
                    if all {
                        let ref_base = rec.target_seq().chars().nth(ref_offset - 1).unwrap_or('-');
                        let alt = Alt {
                            len: 1,
                            vartype: VarType::Null,
                            alt_base: ref_base.to_string(),
                        };
                        let pileup = Pileup {
                            ref_name: ref_name.to_string(),
                            ref_pos,
                            ref_base,
                            alt,
                            uid: uid.clone(),
                        };
                        pileup_vec.push(pileup);
                    }
                }
            }
            'I' => {
                let ref_base = rec.target_seq().chars().nth(ref_offset - 1).unwrap_or('-');
                if ref_base == '-' {
                    ref_offset += len as usize;
                    continue;
                }
                let alt_base = g.iter().map(|(_, c2)| c2).collect::<String>();
                let alt = Alt {
                    len,
                    vartype: VarType::Ins,
                    alt_base,
                };
                let pileup = Pileup {
                    ref_name: ref_name.to_string(),
                    ref_pos,
                    ref_base,
                    alt,
                    uid: uid.clone(),
                };
                pileup_vec.push(pileup);
                ref_offset += len as usize;
            }
            'D' => {
                let ref_base = rec.target_seq().chars().nth(ref_offset - 1).unwrap_or('-');
                if ref_base == '-' {
                    ref_offset += len as usize;
                    ref_pos += len;
                    continue;
                }
                let alt_base = g.iter().map(|(c1, _)| c1).collect::<String>();
                let alt = Alt {
                    len,
                    vartype: VarType::Del,
                    alt_base,
                };
                let pileup = Pileup {
                    ref_name: ref_name.to_string(),
                    ref_pos,
                    ref_base,
                    alt,
                    uid: uid.clone(),
                };
                pileup_vec.push(pileup);
                ref_offset += len as usize;
                ref_pos += len;
            }
            'X' => {
                for i in 0..len {
                    let ref_base = rec.target_seq().chars().nth(ref_offset).unwrap_or('-');
                    let alt_base = g.iter().map(|(_, c2)| c2).nth(i as usize).unwrap_or(&'-');
                    ref_pos += 1;
                    ref_offset += 1;
                    let alt = Alt {
                        len: 1,
                        vartype: VarType::Snp,
                        alt_base: alt_base.to_string(),
                    };
                    let pileup = Pileup {
                        ref_name: ref_name.to_string(),
                        ref_pos,
                        ref_base,
                        alt,
                        uid: uid.clone(),
                    };
                    pileup_vec.push(pileup);
                }
            }
            _ => {}
        }
    }

    Ok(pileup_vec)
}
