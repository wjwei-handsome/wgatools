use std::{
    collections::HashMap,
    fs::File,
    io::{Seek, Write},
};

use itertools::enumerate;
use serde::{Deserialize, Serialize};

use crate::parser::maf::MAFReader;

pub fn build_index(mafreader: &mut MAFReader<File>, idx_wtr: Box<dyn Write>) {
    // init a MAfIndex2 struct
    let mut idx: MafIndex = HashMap::new();

    loop {
        let offset = mafreader.inner.stream_position().unwrap();
        let record = mafreader.records().next();
        if record.is_none() {
            break;
        }
        let record = record.unwrap().unwrap();
        for (ord, sline) in enumerate(record.slines) {
            let name = sline.name;
            let start = sline.start;
            let end = sline.start + sline.align_size;
            let size = sline.size;

            idx.entry(name.clone()).or_insert(MafIndexItem {
                ivls: Vec::new(),
                size,
                ord,
            });
            idx.get_mut(&name)
                .unwrap()
                .ivls
                .push(IvP { start, end, offset });
        }
    }
    serde_json::to_writer(idx_wtr, &idx).unwrap();
}

pub type MafIndex = HashMap<String, MafIndexItem>;

#[derive(Debug, Serialize, Deserialize)]
pub struct MafIndexItem {
    pub ivls: Vec<IvP>,
    size: u64,
    pub ord: usize,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct IvP {
    pub start: u64,
    pub end: u64,
    pub offset: u64,
}
