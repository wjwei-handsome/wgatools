use crate::parser::common::AlignRecord;
use crate::parser::maf::MAFReader;
use noodles_bam as bam;
use noodles_sam as sam;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io;
use std::num::NonZeroUsize;

use noodles_sam::header::record::value::map;
use noodles_sam::header::record::value::map::header::SortOrder;
use noodles_sam::header::record::value::map::Program;
use noodles_sam::header::record::value::{map::ReferenceSequence, Map};

/// Convert a MAF Reader to output a BAM file
pub fn maf2bam<R: io::Read>(mafreader: &mut MAFReader<R>, outputpath: &str) {
    // init a writer with output path, not support stdout
    let mut writer = bam::writer::Builder::default()
        .build_from_path(outputpath)
        .unwrap(); // TODO: handle IO error

    // init a header
    let mut hd = Map::<map::Header>::default();
    *hd.sort_order_mut() = Some(SortOrder::Unsorted);
    let mut header = sam::Header::builder().set_header(hd).build();

    // add @PG header
    let program = Map::<Program>::try_from(vec![
        (String::from("VN"), String::from("1.6.0")),
        (String::from("CL"), String::from("wgalib::maf2bam")),
    ])
    .unwrap();
    header.programs_mut().insert("maf2bam".to_string(), program);

    // collect all maf records and sort by target name and target pos
    let mut mafrecords = mafreader
        .records()
        .map(|rec| rec.unwrap())
        .collect::<Vec<_>>();
    mafrecords.sort();

    // collect all target names and sizes in parallel
    let name_size = mafrecords
        .par_iter()
        .map(|mafrec| (mafrec.target_name(), mafrec.target_length()))
        .collect::<Vec<_>>();
    // define target-id map
    let mut name_id_map = HashMap::<&str, u64>::default();
    let mut id = 0;
    for (name, size) in name_size {
        if !header.reference_sequences().contains_key(name) {
            name_id_map.insert(name, id);
            header.reference_sequences_mut().insert(
                name.parse().unwrap(),
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(size as usize).unwrap()),
            );
            id += 1;
        } else {
            continue;
        }
    }
    // write header
    writer.write_header(&header).unwrap();

    // convert maf records to bam records in parallel
    let bamrecs = mafrecords
        .par_iter()
        .map(|record| record.convert2bam(&name_id_map))
        .collect::<Vec<_>>();

    // write bam records
    bamrecs.into_iter().for_each(|record| {
        writer.write_record(&header, &record).unwrap() // TODO: handle io error
    });
}
