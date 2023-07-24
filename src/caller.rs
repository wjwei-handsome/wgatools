use crate::parser::cigar::cigar_cat_ext;
use crate::parser::common::AlignRecord;
use crate::parser::maf::{MAFReader, MAFRecord};
use itertools::Itertools;
use noodles_vcf::header::record::value::map::{Contig, Info};
use noodles_vcf::record::genotypes::sample::Value;
use noodles_vcf::record::info::field;
use noodles_vcf::record::info::field::key;
use noodles_vcf::record::Genotypes;
use noodles_vcf::{self as vcf, header::record::value::Map, record::Position, Record};
use rayon::prelude::*;
use std::io::BufWriter;

///
///A example:
///
/// ACGATGCTAGCT---ACG
/// AC--TGCTAACTGGGACG
/// within alignment: snp | ins | del | tandem expansion | tandem contraction | Repeat expansion | Repeat contraction
/// between alignment: INS | DEL | Repeat expansion | Repeat contraction

pub fn test_call() {
    let mut reader =
        MAFReader::from_path("/Users/wjwei/NAM2B73v5.Zm-CML333__TO__Zm-B73v5.maf").unwrap();
    let _header = &reader.header;

    let mut wrt = BufWriter::new(std::fs::File::create("test.vcf").unwrap());
    let mut vcf_wrt = vcf::Writer::new(&mut wrt);

    let id = key::SV_LENGTHS;
    let info = Map::<Info>::from(&id);

    let mut header = vcf::Header::builder()
        .add_info(id, info)
        .add_sample_name("sample1")
        .build();
    // vcf_wrt.write_header(&header).expect("TODO: panic message");

    let mut mafrecords = reader
        .records()
        .par_bridge()
        .into_par_iter()
        .map(|rec| rec.unwrap())
        // .for_each(|mut rec| {
        //     println!("got a record");
        // });
        .collect::<Vec<_>>();
    mafrecords.sort();
    let final_var_recs = mafrecords
        .par_iter()
        .map(call_maf_snp_indel_within)
        .flatten()
        .collect::<Vec<_>>();
    let name_size = mafrecords
        .par_iter()
        .map(|mafrec| (mafrec.target_name(), mafrec.target_length()))
        .collect::<Vec<_>>();
    for (name, size) in name_size {
        if !header.contigs().contains_key(name) {
            let mut contigmap = Map::<Contig>::new();
            *contigmap.length_mut() = Some(size as usize);
            header
                .contigs_mut()
                .insert(name.parse().unwrap(), contigmap);
        } else {
            continue;
        }
    }
    // println!("var_rec collected");
    // let header = header_builder.build();
    vcf_wrt.write_header(&header).expect("TODO: panic message");
    let mut snp_count = 0;
    let mut ins_count = 0;
    let mut del_count = 0;
    // let _ = final_var_recs.par_iter()
    //     .map(|varrec| vcf_wrt.write_record(&header, varrec).expect("panics"));
    for mut rec in final_var_recs {
        match rec.info().get(&key::SV_TYPE) {
            Some(value) => {
                if value == Some(&field::Value::String("INS".into())) {
                    let new_id = format!("INS{}", ins_count);
                    *rec.ids_mut() = new_id.parse().unwrap();
                    ins_count += 1;
                } else if value == Some(&field::Value::String("DEL".into())) {
                    let new_id = format!("DEL{}", del_count);
                    *rec.ids_mut() = new_id.parse().unwrap();
                    del_count += 1;
                }
            }
            None => {
                let new_id = format!("SNP{}", snp_count);
                *rec.ids_mut() = new_id.parse().unwrap();
                snp_count += 1;
            }
        }

        vcf_wrt
            .write_record(&header, &rec)
            .expect("TODO: panic message");
    }
}

fn get_variant_rec(
    chro: &str,
    pos: usize,
    ref_base: &str,
    alt_base: &str,
    // id: &str,
    info: &str,
) -> Record {
    let mut record = Record::builder()
        .set_chromosome(chro.parse().unwrap())
        .set_position(Position::from(pos))
        .set_reference_bases(ref_base.parse().unwrap())
        .set_alternate_bases(alt_base.parse().unwrap())
        .set_info(info.parse().unwrap())
        // .set_ids(id.parse().unwrap())
        .build()
        .unwrap();
    let keys = "GT".parse().unwrap();
    let genotypes = Genotypes::new(keys, vec![vec![Some(Value::from("1|1"))]]);

    *record.genotypes_mut() = genotypes;
    record
}

fn call_maf_snp_indel_within(mafrec: &MAFRecord) -> Vec<Record> {
    // target:ACG-TTTGATGCTAGCT---ACG
    // query :ACCATTT--TGCTAACTGGGACG

    let mut var_recs = Vec::new();

    let mut target_current_offset = mafrec.target_start();
    let mut query_current_offset = mafrec.query_start();

    let chro = mafrec.target_name();

    let t_start = mafrec.target_start();
    let q_start = mafrec.query_start();

    let mut t_seq_ref = mafrec.target_seq().to_string();
    t_seq_ref.retain(|c| c != '-');

    let mut q_seq_ref = mafrec.query_seq().to_string();
    q_seq_ref.retain(|c| c != '-');

    let t_seq_iter = mafrec.target_seq().chars();
    let q_seq_iter = mafrec.query_seq().chars();
    let group_by_iter = t_seq_iter
        .zip(q_seq_iter)
        .group_by(|(c1, c2)| cigar_cat_ext(c1, c2));

    for (k, g) in group_by_iter.into_iter() {
        let len = g.count() as u64;
        match k {
            '=' => {
                target_current_offset += len;
                query_current_offset += len;
            }
            'I' => {
                let t_slice_start = (target_current_offset - t_start - 1) as usize;
                let t_slice_end = t_slice_start + 1;

                let q_slice_start = (query_current_offset - q_start - 1) as usize;
                let q_slice_end = q_slice_start + len as usize + 1;

                let info = format!("SVTYPE=INS;SVLEN={}", len);
                // let id = format!("INS{}", ins_count);
                let ref_base = &t_seq_ref[t_slice_start..t_slice_end];
                let alt_base = &q_seq_ref[q_slice_start..q_slice_end];
                let record = get_variant_rec(
                    chro,
                    target_current_offset as usize,
                    ref_base,
                    alt_base,
                    // &id,
                    &info,
                );

                // vcf_wrt.write_record(&header, &record).expect("TODO: panic message");
                var_recs.push(record);
                query_current_offset += len;
            }
            'D' => {
                let t_slice_start = (target_current_offset - t_start - 1) as usize;
                let t_slice_end = t_slice_start + len as usize + 1;

                let q_slice_start = (query_current_offset - q_start - 1) as usize;
                let q_slice_end = q_slice_start + 1;

                let info = format!("SVTYPE=DEL;SVLEN={}", len);
                // let id = format!("DEL{}", del_count);
                let ref_base = &t_seq_ref[t_slice_start..t_slice_end];
                let alt_base = &q_seq_ref[q_slice_start..q_slice_end];
                let record = get_variant_rec(
                    chro,
                    target_current_offset as usize,
                    ref_base,
                    alt_base,
                    // &id,
                    &info,
                );
                // vcf_wrt.write_record(&header, &record).expect("TODO: panic message");
                var_recs.push(record);
                target_current_offset += len;
            }
            'X' => {
                // println!("{} {} {}", current_offset, len, k);
                let t_slice_start = (target_current_offset - t_start) as usize;
                let t_slice_end = t_slice_start + 1;

                let q_slice_start = (query_current_offset - q_start) as usize;
                let q_slice_end = q_slice_start + 1;

                let info = "TYPE=SNP".to_string();
                // let id = format!("SNP{}", snp_count);
                let ref_base = &t_seq_ref[t_slice_start..t_slice_end];
                let alt_base = &q_seq_ref[q_slice_start..q_slice_end];
                let record = get_variant_rec(
                    chro,
                    target_current_offset as usize,
                    ref_base,
                    alt_base,
                    // &id,
                    &info,
                );
                // vcf_wrt.write_record(&header, &record).expect("TODO: panic message");
                var_recs.push(record);
                // println!("{}\t{}\t{}\t{}\tSNP\t{}\t{}",
                //          target_current_offset+1,
                //          target_current_offset+len+1,
                //          query_current_offset+1,
                //          query_current_offset+len+1,
                //          &t_seq_ref[t_slice_start..t_slice_end],
                //          &q_seq_ref[q_slice_start..q_slice_end],
                //          );
                query_current_offset += len;
                target_current_offset += len;
            }
            _ => {}
        }
    }
    println!("done a maf record");
    var_recs
}
