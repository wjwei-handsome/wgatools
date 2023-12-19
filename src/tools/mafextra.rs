use regex::Regex;
use serde::{Deserialize, Serialize};
use std::cmp::{max, min};
use std::fmt::Display;
use std::fs::File;
use std::io::{BufReader, Write};

use super::index::{IvP, MafIndex};

use crate::errors::{ParseGenomeRegionErrKind, WGAError};
use crate::parser::maf::{MAFReader, MAFWriter};
use crate::utils::parse_str2u64;
use csv::ReaderBuilder;
use std::io::Read;
use std::io::Seek;

use rust_lapper::{Interval, Lapper};

// fn maf_extract_iter<R: Read>(
//     _regions: &Option<Vec<String>>,
//     _region_file: &Option<String>,
//     _outputpath: &str,
//     _mafreader: &mut MAFReader<R>,
// ) {
//     todo!()
// }

pub fn maf_extract_idx<R: Read + Send + Seek>(
    regions: &Option<Vec<String>>,
    region_file: &Option<String>,
    mafreader: &mut MAFReader<R>,
    mafindex: MafIndex,
    writer: &mut dyn Write,
) -> Result<Vec<GenomeRegion>, WGAError> {
    let input_regions = get_input_regions(regions, region_file)?;
    let mut sub_maf_wtr = MAFWriter::new(writer);
    let header = "#maf version=1.6 cmd=maf_extract";
    sub_maf_wtr.write_header(header.to_owned())?;
    let failed_regions =
        extract_sub_blocks_with_idx(mafindex, input_regions, mafreader, &mut sub_maf_wtr)?;
    Ok(failed_regions)
}

fn get_input_regions(
    regions: &Option<Vec<String>>,
    region_file: &Option<String>,
) -> Result<Vec<GenomeRegion>, WGAError> {
    // judge regions and region_file
    // acutally it's unnecessary
    if regions.is_none() && region_file.is_none() {
        return Err(WGAError::EmptyRegion);
    }

    // init input regions
    let mut input_regions = Vec::new();

    // read input region_vec
    if let Some(regions) = regions {
        for region in regions {
            let genome_region = GenomeRegion::try_from(region.to_string())?;
            input_regions.push(genome_region);
        }
    }

    // read input region_file
    if let Some(region_file) = region_file {
        let reader = BufReader::new(File::open(region_file)?);
        let regions = read_genome_region(reader)?;
        input_regions.extend(regions);
    }
    Ok(input_regions)
}

#[derive(Debug, Serialize, Deserialize)]
pub struct GenomeRegion {
    name: String,
    start: u64,
    end: u64,
}

impl TryFrom<String> for GenomeRegion {
    type Error = WGAError;
    fn try_from(value: String) -> Result<Self, Self::Error> {
        let re = Regex::new(r"^([a-zA-Z0-9.@_-]+):([0-9]+)-([0-9]+)$")?;
        match re.captures(&value) {
            Some(caps) => {
                let name = caps
                    .get(1)
                    .ok_or(WGAError::UnexceptedRegexError("name".to_string()))?
                    .as_str()
                    .to_string();
                let start = caps
                    .get(2)
                    .ok_or(WGAError::UnexceptedRegexError("start".to_string()))?
                    .as_str();
                let start = parse_str2u64(start)?;
                let end = caps
                    .get(3)
                    .ok_or(WGAError::UnexceptedRegexError("end".to_string()))?
                    .as_str();
                let end = parse_str2u64(end)?;
                if start > end {
                    return Err(WGAError::ParseGenomeRegion(
                        ParseGenomeRegionErrKind::StartGTEnd(start, end),
                    ));
                }
                Ok(GenomeRegion { name, start, end })
            }
            None => Err(WGAError::ParseGenomeRegion(
                ParseGenomeRegionErrKind::FormatNotMatch(value),
            )),
        }
    }
}

impl Display for GenomeRegion {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}:{}-{}", self.name, self.start, self.end)
    }
}

fn read_genome_region<R: Read>(reader: R) -> Result<Vec<GenomeRegion>, WGAError> {
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(reader);
    let mut regions = Vec::new();
    for result in rdr.deserialize() {
        let record: GenomeRegion = result?;
        if record.start > record.end {
            return Err(WGAError::ParseGenomeRegion(
                ParseGenomeRegionErrKind::StartGTEnd(record.start, record.end),
            ));
        }
        regions.push(record);
    }
    Ok(regions)
}

type Iv = Interval<u64, u64>;

fn ivp2iv(ivp: &IvP) -> Iv {
    Iv {
        start: ivp.start,
        stop: ivp.end,
        val: ivp.offset,
    }
}

fn extract_sub_blocks_with_idx<R: Read + Send + Seek, W: Write>(
    mafidx: MafIndex,
    regions: Vec<GenomeRegion>,
    mafreader: &mut MAFReader<R>,
    mafwriter: &mut MAFWriter<W>,
) -> Result<Vec<GenomeRegion>, WGAError> {
    let mut failed_regions = Vec::new();
    // TODO: parallel genearte sub-maf-blocks
    for givl in regions.into_iter() {
        match mafidx.get(&givl.name) {
            Some(item) => {
                let hit_ivps = &item.ivls;
                let hit_givls = hit_ivps.iter().map(ivp2iv).collect::<Vec<Iv>>();
                let lapper = Lapper::new(hit_givls);
                let find = lapper.find(givl.start, givl.end).collect::<Vec<&Iv>>();
                let find_num = find.len();
                let ord = item.ord;
                match find_num {
                    0 => {
                        failed_regions.push(givl);
                        continue;
                    }
                    _ => {
                        for block in find {
                            let offset = block.val;
                            mafreader.inner.seek(std::io::SeekFrom::Start(offset))?;
                            let mut mafrec =
                                mafreader.records().next().ok_or(WGAError::EmptyRecord)??;

                            let b_start = block.start;
                            let b_end = block.stop;

                            let g_start = givl.start;
                            let g_end = givl.end;

                            if g_start <= b_start && g_end >= b_end {
                                mafwriter.write_record(&mafrec)?;
                                continue;
                            }

                            let r_start = max(b_start, g_start);
                            let r_end = min(b_end, g_end);

                            mafrec.slice_block(r_start, r_end, ord);

                            mafwriter.write_record(&mafrec)?;
                        }
                    }
                }
            }
            None => {
                failed_regions.push(givl);
                continue;
            }
        };
    }
    Ok(failed_regions)
}
