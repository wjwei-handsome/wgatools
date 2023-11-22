use crate::parser::maf::MAFReader;
use crate::utils::output_writer;
use noodles::sam::header::record::value::map;
use noodles::sam::header::record::value::map::header::SortOrder;
use noodles::sam::record::ReadName;
use noodles::sam::{
    self as sam,
    header::record::value::{
        map::{Program, ReferenceSequence},
        Map,
    },
};

use std::io;
use std::num::NonZeroUsize;

/// Convert a MAF Reader to output a BAM file
// pub fn maf2bam<R: io::Read + Send>(mafreader: &mut MAFReader<R>, outputpath: &str) {
//     // init a writer with output path, not support stdout
//     let mut writer = bam::writer::Builder::default()
//         .build_from_path(outputpath)
//         .unwrap(); // TODO: handle IO error
//
//     // init a header
//     let mut hd = Map::<map::Header>::default();
//     *hd.sort_order_mut() = Some(SortOrder::Unsorted);
//     let mut header = sam::Header::builder().set_header(hd).build();
//
//     // add @PG header
//     let program = Map::<Program>::try_from(vec![
//         (String::from("VN"), String::from("1.6.0")),
//         (String::from("CL"), String::from("wgalib::maf2bam")),
//     ])
//     .unwrap();
//     header.programs_mut().insert("maf2bam".to_string(), program);
//
//     // collect all maf records and sort by target name and target pos
//     let mafrecords = mafreader
//         .records()
//         .par_bridge()
//         .map(|rec| rec.unwrap())
//         .collect::<Vec<_>>();
//     // mafrecords.sort();
//
//     // collect all target names and sizes in parallel
//     let name_size = mafrecords
//         .par_iter()
//         .map(|mafrec| (mafrec.target_name(), mafrec.target_length()))
//         .collect::<Vec<_>>();
//     // define target-id map
//     let mut name_id_map = HashMap::<&str, u64>::default();
//     let mut id = 0;
//     for (name, size) in name_size {
//         if !header.reference_sequences().contains_key(name) {
//             name_id_map.insert(name, id);
//             header.reference_sequences_mut().insert(
//                 name.parse().unwrap(),
//                 Map::<ReferenceSequence>::new(NonZeroUsize::try_from(size as usize).unwrap()),
//             );
//             id += 1;
//         } else {
//             continue;
//         }
//     }
//     // write header
//     writer.write_header(&header).unwrap();
//
//     // convert maf records to bam records in parallel
//     let bamrecs = mafrecords
//         .par_iter()
//         .map(|record| record.convert2bam(&name_id_map))
//         .collect::<Vec<_>>();
//
//     // write bam records
//     bamrecs.into_iter().for_each(|record| {
//         match writer.write_record(&header, &record) {
//             Ok(_) => {println!("{:?}",record.cigar())}
//             Err(e) => {
//                 println!("{}:{:?}", e,record.cigar());
//             }
//         } // TODO: handle io error
//     });
// }

pub fn maf2sam<R: io::Read + Send>(_mafreader: &mut MAFReader<R>, outputpath: &str) {
    let writer = output_writer(outputpath);
    let mut sam_writer = sam::Writer::new(writer);
    let mut header = Map::<map::Header>::default();
    *header.sort_order_mut() = Some(SortOrder::Unsorted);
    let header = sam::Header::builder()
        .set_header(header)
        .add_reference_sequence(
            "sq0".parse().unwrap(),
            Map::<ReferenceSequence>::new(NonZeroUsize::try_from(8).unwrap()),
        )
        .add_reference_sequence(
            "sq1".parse().unwrap(),
            Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13).unwrap()),
        )
        .add_reference_sequence(
            "sq2".parse().unwrap(),
            Map::<ReferenceSequence>::new(NonZeroUsize::try_from(21).unwrap()),
        )
        .add_program("noodles-sam", Map::<Program>::default())
        .add_comment("an example SAM written by noodles-sam")
        .build();

    sam_writer.write_header(&header).unwrap();
    let mut record = sam::alignment::Record::default();
    let read_name: ReadName = "sq2".parse().unwrap();
    *record.read_name_mut() = Some(read_name);
    sam_writer.write_record(&header, &record).unwrap();

    // init a writer with output path, not support stdout
    // let mut bamheader = rust_htslib::bam::header::Header::new();

    // collect all maf records and sort by target name and target pos
    // let mut mafrecords = mafreader
    //     .records()
    //     .par_bridge()
    //     .map(|rec| rec.unwrap())
    //     .collect::<Vec<_>>();
    // mafrecords.sort();

    // collect all target names and sizes in parallel
    // let name_size = mafrecords
    //     .par_iter()
    //     .map(|mafrec| (mafrec.target_name(), mafrec.target_length()))
    //     .collect::<Vec<_>>();
    // let mut name_id_map = HashMap::<&str, u64>::default();
    // let mut id = 0;
    // let mut bamheader_sample_rec = rust_htslib::bam::header::HeaderRecord::new(b"SQ");
    // let bamheader_sample_rec = bamheader_sample_rec
    //     .push_tag(b"SN", "chr1")
    //     .push_tag(b"LN", 1000);
    // bamheader.push_record(&bamheader_sample_rec);
    // // for (name, size) in name_size {
    //     if name_id_map.contains_key(name) {
    //         continue;
    //     } else {
    //         let bamheader_sample_rec = bamheader_sample_rec.push_tag(b"SN", name)
    //             .push_tag(b"LN", size);
    //         bamheader.push_record(&bamheader_sample_rec);
    //         name_id_map.insert(name, id);
    //         id += 1;
    //     }
    // }
    // let mut writer = BamWriter::from_path(outputpath,&bamheader, Format::Sam).unwrap(); // TODO: handle IO error

    // convert maf records to bam records in parallel
    // let bamrecs = mafrecords
    //     .par_iter()
    //     .map(|record| convert2bam2(record, &name_id_map))
    //     .collect::<Vec<_>>();

    // write bam records
    // bamrecs.into_iter().for_each(|record| {
    //     writer.write(&record).unwrap(); // TODO: handle io error
    // });
    // let mut std_wtr = BamWriter::from_stdout(&bamheader, Format::Sam).unwrap();
    // let mut bamrec = BamRecord::default();
    // bamrec.set_pos(1);
    // bamrec.set_tid(0);
    // bamrec.set_flags(0);
    // let seq = b"ATGC";
    // let quals = vec![255 as u8; 4];
    // let fake_cigar_vec = vec![Cigar::HardClip(289573083)];
    // let fake_cigar = CigarString(fake_cigar_vec);
    // bamrec.set(b"qname", Some(&fake_cigar), seq, &quals);
    // bamrec.set_mapq(255);
    // bamrec.set_mpos(0);
    // bamrec.set_mtid(-1);
    // std_wtr.write(&bamrec).unwrap();
    // println!("{:?}", bamrec);
    // match writer.write(&bamrec) {
    //     Ok(_) => {}
    //     Err(e) => {
    //         println!("{}", e);
    //     }
    // }
}

// fn convert2bam2(maferc: &MAFRecord, name_id_map: &HashMap<&str, u64>) -> BamRecord {
//     // init a bam record
//     let mut bamrec = BamRecord::new();

//     // set bam record query name
//     let q_name = maferc.query_name();
//     // bamrec.set_qname(q_name.as_bytes());

//     // set bam record flags: it always in empty in whole genome alignment
//     bamrec.set_flags(0);

//     bamrec.set_mapq(255);
//     bamrec.set_mpos(0);
//     bamrec.set_mtid(-1);
//     // set bam record reference sequence id ref to name-id-map
//     let t_name = maferc.target_name();
//     let t_id = *name_id_map.get(t_name).unwrap();
//     bamrec.set_tid(t_id as i32);

//     // set bam record position
//     let t_start = maferc.target_start(); // 0-based to 1-based
//     bamrec.set_pos(t_start as i64);

//     // set bam record cigar
//     let gen_cigar = parse_maf_seq_to_cigar(maferc, true);
//     let fuckcigar_vec = vec![Cigar::HardClip(289573083)];
//     let cigar = CigarString(fuckcigar_vec);
//     // let cigar = gen_cigar.bamcigar;
//     println!("cigar0{:?}", cigar[0]);

//     // set bam record sequence
//     let mut q_seq_ref = maferc.query_seq().to_string();
//     q_seq_ref.retain(|c| c != '-'); // remove gap;should be UPPER?

//     // let q_seq_ref = b"A";

//     let qual = vec![255 as u8; q_seq_ref.len()];

//     bamrec.set(q_name.as_bytes(), Some(&cigar), q_seq_ref.as_bytes(), &qual);

//     // set bam record tags
//     let edit_dist = gen_cigar.mismatch_count + gen_cigar.ins_count + gen_cigar.del_count;
//     let nm_tag = String::from("NM:i:") + &*edit_dist.to_string();

//     // return bam record
//     bamrec
// }
