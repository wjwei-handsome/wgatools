use crate::{
    errors::WGAError,
    parser::{
        common::recount_align_size,
        maf::{MAFReader, MAFRecord, MAFSLine, MAFWriter},
    },
};
use std::io::{Read, Write};

// mian function of chunk maf
// A 0 5  + 5  -----ATCGT
// B 0 10 + 10 GGGGGATCGT

// chunk size = 5

// A 0 0 + 0 -----
// B 0 5 + 0 GGGGG

// A 0 5 + 5 ATCGT
// B 5 5 + 5 ATCGT
pub fn chunk_maf<R: Read + Send>(
    mut reader: MAFReader<R>,
    chunk_length: u64,
    writer: &mut dyn Write,
) -> Result<(), WGAError> {
    // init a MAFWriter
    let mut mafwtr = MAFWriter::new(writer);
    // write header
    let header = format!("#maf version=1.6 split_length={}", chunk_length);
    mafwtr.write_header(header)?;

    // chunk each block
    for rec in reader.records() {
        let rec = rec?;
        let block_length = rec.slines[0].seq.len() as u64;

        // init sline_end_vec
        let mut sline_end_vec = vec![];
        for s in &rec.slines {
            sline_end_vec.push(s.start);
        }
        // chunk each s-line in a block until the end
        let mut chunk_start = 0;
        let mut chunk_end = chunk_length;

        // chunk each s-line in a block until the end
        while chunk_end < block_length {
            let new_rec = emit_new_maf_rec(&rec, chunk_start, chunk_end, &mut sline_end_vec)?;
            mafwtr.write_record(&new_rec)?;
            chunk_start = chunk_end;
            chunk_end += chunk_length;
        }

        // last chunk
        let new_rec = emit_new_maf_rec(&rec, chunk_start, block_length, &mut sline_end_vec)?;
        mafwtr.write_record(&new_rec)?;
    }

    Ok(())
}

// emit new maf rec
fn emit_new_maf_rec(
    rec: &MAFRecord,
    chunk_start: u64,
    chunk_end: u64,
    end_vec: &mut [u64],
) -> Result<MAFRecord, WGAError> {
    let mut new_rec = MAFRecord {
        score: rec.score,
        slines: vec![],
        query_idx: 1,
    };
    for (i, sline) in rec.slines.iter().enumerate() {
        let new_seq = &sline.seq[chunk_start as usize..chunk_end as usize];
        let (align_size, _) = recount_align_size(new_seq);
        let new_sline = MAFSLine {
            mode: 's',
            name: sline.name.clone(),
            start: end_vec[i],
            align_size,
            strand: sline.strand,
            size: sline.size,
            seq: new_seq.to_string(),
        };
        new_rec.slines.push(new_sline);
        end_vec[i] += align_size;
    }
    Ok(new_rec)
}
