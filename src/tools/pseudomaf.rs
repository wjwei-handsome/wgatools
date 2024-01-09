use crate::{
    errors::WGAError,
    parser::{
        cigar::gen_pesudo_maf_by_cigar,
        common::AlignRecord,
        paf::{PAFReader, PafRecord},
    },
};
use std::{
    collections::HashMap,
    io::{Read, Write},
};

// main function of generate pesudo MAF from PAF

pub fn generate_pesudo_maf<R: Read + Send>(
    mut reader: PAFReader<R>,
    writer: &mut dyn Write,
) -> Result<(), WGAError> {
    writeln!(writer, "##maf version=1")?;

    //  groupby the target name and query name and sort by target start
    let mut target_qeury_recs_map: HashMap<String, HashMap<String, Vec<PafRecord>>> =
        HashMap::new();
    for pafrec in reader.records() {
        let rec = pafrec?;
        let target_name = rec.target_name().to_string();
        let query_name = rec.query_name().to_string();
        let target_start = rec.target_start();
        let query_recs_map = target_qeury_recs_map
            .entry(target_name)
            .or_insert(HashMap::new());
        let query_rec_vec = query_recs_map.entry(query_name).or_insert(Vec::new());
        // sort by target start
        let idx = query_rec_vec
            .binary_search_by(|probe| probe.target_start().cmp(&target_start))
            .unwrap_or_else(|e| e);
        query_rec_vec.insert(idx, rec);
    }

    // 2. merge the rec_vec
    let final_map = target_qeury_recs_map
        .into_iter()
        .map(|(target_name, query_recs_map)| {
            let query_recs_vec = merge_query_recs_map(query_recs_map)?;
            Ok::<_, WGAError>((target_name, query_recs_vec))
        })
        .collect::<Result<HashMap<_, _>, WGAError>>()?;

    // 2. for each target, generate a pesudo maf
    for (target_name, rec_vec) in final_map {
        // println!("{:?}", rec_vec);
        write_pmaf(writer, rec_vec, &target_name)?;
    }

    Ok(())
}

// merge sorted rec_vec to a single rec
fn merge_query_recs_map(
    query_recs_map: HashMap<String, Vec<PafRecord>>,
) -> Result<Vec<PafRecord>, WGAError> {
    let mut pafrec_vec = Vec::new();
    for (_, rec_vec) in query_recs_map {
        if rec_vec.len() == 1 {
            pafrec_vec.extend(rec_vec);
        } else {
            // merge the rec_vec to a single rec
            // new target_start is the min of target_start
            // new target_end is the max of target_end
            // new cigar is the merge of cigar, padding by `{n}D`, n is the distance between two rec

            // take first element as the new rec
            let mut rec_vec_iter = rec_vec.into_iter();
            let mut new_rec = match rec_vec_iter.next() {
                Some(rec) => rec,
                None => return Err(WGAError::Other(anyhow::anyhow!("rec_vec is empty"))),
            };
            let mut new_cigar = String::new();
            new_cigar.push_str(new_rec.get_cigar_str()?);
            let mut last_target_end = new_rec.target_end();
            for rec in rec_vec_iter {
                let target_start = rec.target_start();
                let target_end = rec.target_end();
                // trim "cg:Z:"
                let pure_cigar = rec.get_cigar_str()?.trim_start_matches("cg:Z:");
                let gap_len = target_start - last_target_end;
                new_cigar.push_str(&format!("{}D", gap_len));
                new_cigar.push_str(pure_cigar);
                last_target_end = target_end;
            }

            new_rec.tags = vec![new_cigar];
            new_rec.target_end = last_target_end;
            pafrec_vec.push(new_rec);
        }
    }
    Ok(pafrec_vec)
}

fn write_pmaf(
    writer: &mut dyn Write,
    rec_vec: Vec<PafRecord>,
    target_name: &str,
) -> Result<(), WGAError> {
    writeln!(writer, "a score=0")?;
    let mut first_flag = true;

    for rec in rec_vec {
        let target_size = rec.target_length();
        let target_start = rec.target_start();
        if first_flag {
            write!(
                writer,
                "s\t{}\t0\t{}\t+\t{}\t",
                target_name, target_size, target_size
            )?;
            for _ in 0..target_size {
                write!(writer, "N")?;
            }
            writeln!(writer)?;
        }
        first_flag = false;
        let query_name = rec.query_name();
        let query_size = rec.query_length();
        let cigar = rec.get_cigar_str()?;

        write!(
            writer,
            "s\t{}\t0\t{}\t+\t{}\t",
            query_name, query_size, query_size
        )?;
        for _ in 0..target_start {
            write!(writer, "-")?;
        }
        let gen_seq_line = gen_pesudo_maf_by_cigar(cigar)?;
        write!(writer, "{}", gen_seq_line)?;
        let seq_line_len = gen_seq_line.len() as u64;
        for _ in 0..(target_size - seq_line_len - target_start) {
            write!(writer, "-")?;
        }
        writeln!(writer)?;
    }
    writeln!(writer)?;
    Ok(())
}
