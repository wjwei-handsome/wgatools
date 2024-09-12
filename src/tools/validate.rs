use crate::{
    errors::WGAError,
    parser::{
        common::AlignRecord,
        paf::{PAFReader, PafRecord},
    },
};
use rayon::prelude::*;
use std::fmt;
use std::io::{Read, Write};

/// Check query&target start&end position by CIGAR
/// query_start + Match/Mismatch + INS_size = query_end
/// ref_start + Match/Mismatch + DEL_size = ref_end

#[derive(Default)]
struct Validations {
    total: usize,
    query_invalid: usize,
    query_inv_list: Vec<String>,
    ref_invalid: usize,
    ref_inv_list: Vec<String>,
    fix_paf_recs: Vec<PafRecord>,
}

impl fmt::Display for Validations {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Total records: {}", self.total)?;
        writeln!(f, "Query invalid records: {}", self.query_invalid)?;
        writeln!(f, "Target invalid records: {}", self.ref_invalid)?;
        writeln!(f, "Query invalid list:")?;
        for query_uid in self.query_inv_list.iter() {
            writeln!(f, "{}", query_uid)?;
        }
        writeln!(f, "Target invalid list:")?;
        for ref_uid in self.ref_inv_list.iter() {
            writeln!(f, "{}", ref_uid)?;
        }
        Ok(())
    }
}

/// Validate PAF records in parallel
pub fn parallel_validatepaf<R: Read + Send>(
    mut reader: PAFReader<R>,
    writer: &mut dyn Write,
    fix_writer: Option<Box<dyn Write>>,
    fix_flag: bool,
) -> Result<(), WGAError> {
    let validations = reader
        .records()
        .par_bridge()
        .try_fold(Validations::default, |vd, rec| {
            let rec = rec?;
            process_record(vd, rec, fix_flag)
        })
        .try_reduce(Validations::default, |mut vd1, vd2| {
            vd1.total += vd2.total;
            vd1.query_invalid += vd2.query_invalid;
            vd1.query_inv_list.extend(vd2.query_inv_list);
            vd1.ref_invalid += vd2.ref_invalid;
            vd1.ref_inv_list.extend(vd2.ref_inv_list);
            vd1.fix_paf_recs.extend(vd2.fix_paf_recs);
            Ok(vd1)
        });
    process_validations(validations?, writer, fix_writer)?;
    Ok(())
}

/// process record
fn process_record(
    mut vd: Validations,
    mut rec: PafRecord,
    fix_flag: bool,
) -> Result<Validations, WGAError> {
    vd.total += 1;
    let rec_stat = rec.get_stat().unwrap();

    // check query end
    let exp_query_end = rec.query_start()
        + rec_stat.matched as u64
        + rec_stat.mismatched as u64
        + rec_stat.ins_size as u64;
    if exp_query_end != rec.query_end() {
        vd.query_invalid += 1;
        let query_uid = format!(
            "{}:{}-{}",
            rec.query_name(),
            rec.query_start(),
            rec.query_end()
        );
        vd.query_inv_list.push(query_uid);
        rec.query_end = exp_query_end;
    }

    // check ref end
    let exp_ref_end = rec.target_start()
        + rec_stat.matched as u64
        + rec_stat.mismatched as u64
        + rec_stat.del_size as u64;
    if exp_ref_end != rec.target_end() {
        vd.ref_invalid += 1;
        let ref_uid = format!(
            "{}:{}-{}",
            rec.target_name(),
            rec.target_start(),
            rec.target_end()
        );
        vd.ref_inv_list.push(ref_uid);
        rec.target_end = exp_ref_end;
    }

    if fix_flag {
        vd.fix_paf_recs.push(rec);
    }

    Ok(vd)
}

/// output validations
fn process_validations(
    validations: Validations,
    writer: &mut dyn Write,
    fix_writer: Option<Box<dyn Write>>,
) -> Result<(), WGAError> {
    writeln!(writer, "{}", validations)?;
    // write fix output
    if let Some(writer) = fix_writer {
        let mut pafwtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_writer(writer);
        for rec in validations.fix_paf_recs {
            pafwtr.serialize(rec)?;
        }
    }
    Ok(())
}
