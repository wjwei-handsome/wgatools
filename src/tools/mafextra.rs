use log::error;
use serde::{Deserialize, Serialize};
use std::cmp::{max, min};
use std::fmt::Display;
use std::fs::File;
use std::io::{BufReader, Error, Write};

use super::index::{IvP, MafIndex};
use crate::parser::common::AlignRecord;
use crate::parser::maf::{MAFReader, MAFRecord, MAFWriter};
use csv::{DeserializeRecordsIter, ReaderBuilder};
use rayon::prelude::*;
use std::io::Read;
use std::io::Seek;

use rust_lapper::{Interval, Lapper};

pub fn maf_extract_iter<R: Read>(
    regions: &Option<Vec<String>>,
    region_file: &Option<String>,
    outputpath: &str,
    mafreader: &mut MAFReader<R>,
) {
    todo!()
}

pub fn maf_extract_idx<R: Read + Send + Seek>(
    regions: &Option<Vec<String>>,
    region_file: &Option<String>,
    mafreader: &mut MAFReader<R>,
    mafindex: MafIndex,
    writer: &mut dyn Write,
) {
    let input_regions = get_input_regions(regions, region_file);
    let mut sub_maf_wtr = MAFWriter::new(writer);
    let header = "#maf version=1.6 cmd=maf_extract";
    sub_maf_wtr.write_header(header.to_owned());
    extract_sub_blocks_with_idx(mafindex, input_regions, mafreader, &mut sub_maf_wtr);
}

fn get_input_regions(
    regions: &Option<Vec<String>>,
    region_file: &Option<String>,
) -> Vec<GenomeRegion> {
    // judge regions and region_file
    if regions.is_none() && region_file.is_none() {
        error!("regions or region_file must be specified");
        std::process::exit(1);
    }

    // init input regions
    let mut input_regions = Vec::new();

    // read input region_vec
    if let Some(regions) = regions {
        for region in regions {
            let genome_region = match GenomeRegion::try_from(region.to_string()) {
                Ok(genome_region) => genome_region,
                Err(err) => {
                    error!("region format error: {} in {}", err, region);
                    std::process::exit(1);
                }
            };
            input_regions.push(genome_region);
        }
    }

    // read input region_file
    if let Some(region_file) = region_file {
        let reader = match File::open(region_file) {
            Ok(file) => BufReader::new(file),
            Err(err) => {
                error!("failed to open region file: {}", err);
                std::process::exit(1);
            }
        };
        let regions = read_genome_region(reader);
        match regions {
            Ok(regions) => {
                input_regions.extend(regions);
            }
            Err(err) => {
                error!("region file format error: {}", err);
                std::process::exit(1);
            }
        }
    }
    // println!("{:?}", input_regions);
    input_regions
}

#[derive(Debug, Serialize, Deserialize)]
struct GenomeRegion {
    name: String,
    start: u64,
    end: u64,
}

impl TryFrom<String> for GenomeRegion {
    type Error = Error;
    fn try_from(value: String) -> Result<Self, Self::Error> {
        // TODO: use regex to check the format;
        // TODO: or use trans err.
        let mut iter = value.split(':');
        let name = match iter.next() {
            Some(name) => name.to_string(),
            None => {
                return Err(Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "region format error",
                ))
            }
        };
        let mut iter = match iter.next() {
            Some(start_end) => start_end.split('-'),
            None => {
                return Err(Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "region format error",
                ))
            }
        };
        let start = match iter.next() {
            Some(start) => match start.parse::<u64>() {
                Ok(start) => start,
                Err(_) => {
                    return Err(Error::new(
                        std::io::ErrorKind::InvalidInput,
                        "region format error",
                    ))
                }
            },
            None => {
                return Err(Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "region format error",
                ))
            }
        };
        let end = match iter.next() {
            Some(end) => match end.parse::<u64>() {
                Ok(end) => end,
                Err(_) => {
                    return Err(Error::new(
                        std::io::ErrorKind::InvalidInput,
                        "region format error",
                    ))
                }
            },
            None => {
                return Err(Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "region format error",
                ))
            }
        };
        if start > end {
            return Err(Error::new(
                std::io::ErrorKind::InvalidInput,
                "region format error",
            ));
        }
        Ok(GenomeRegion { name, start, end })
    }
}

impl Display for GenomeRegion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}-{}", self.name, self.start, self.end)
    }
}

fn read_genome_region<R: Read>(reader: R) -> Result<Vec<GenomeRegion>, csv::Error> {
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(reader);
    let mut regions = Vec::new();
    for result in rdr.deserialize() {
        let record: GenomeRegion = result?;
        if record.start > record.end {
            return Err(csv::Error::from(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "region format error",
            ))); // TODO: should be organized
        }
        regions.push(record);
    }
    Ok(regions)
}

fn extract_sub_blocks_with_idx<R: Read + Send + Seek, W: Write>(
    mafidx: MafIndex,
    regions: Vec<GenomeRegion>,
    mafreader: &mut MAFReader<R>,
    mafwriter: &mut MAFWriter<W>,
) {
    let stderr = std::io::stderr();
    for givl in regions.iter() {
        match mafidx.get(&givl.name) {
            Some(item) => {
                let hit_ivps = &item.ivls;
                let hit_givls = hit_ivps.iter().map(IvP2Iv).collect::<Vec<Iv>>();
                let lapper = Lapper::new(hit_givls);
                let find = lapper.find(givl.start, givl.end).collect::<Vec<&Iv>>();
                let find_num = find.len();
                let ord = item.ord;
                match find_num {
                    0 => {
                        stderr
                            .lock()
                            .write_all(format!("region `{}` not found\n", givl,).as_bytes())
                            .unwrap();
                        continue;
                    }
                    _ => {
                        for block in find {
                            let offset = block.val;
                            mafreader
                                .inner
                                .seek(std::io::SeekFrom::Start(offset))
                                .unwrap();
                            let mut mafrec = mafreader.records().next().unwrap().unwrap();

                            let b_start = block.start;
                            let b_end = block.stop;

                            let g_start = givl.start;
                            let g_end = givl.end;

                            if g_start <= b_start && g_end >= b_end {
                                mafwriter.write_record(&mafrec);
                                continue;
                            }

                            let r_start = max(b_start, g_start);
                            let r_end = min(b_end, g_end);

                            mafrec.slice_block(r_start, r_end, ord);

                            mafwriter.write_record(&mafrec);
                        }
                    }
                }
            }
            None => {
                stderr
                    .lock()
                    .write_all(format!("sequence `{}` not found\n", givl.name).as_bytes())
                    .unwrap();
                continue;
            }
        };
    }
}

type Iv = Interval<u64, u64>;

fn IvP2Iv(ivp: &IvP) -> Iv {
    Iv {
        start: ivp.start,
        stop: ivp.end,
        val: ivp.offset,
    }
}
