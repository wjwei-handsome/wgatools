use super::*;
use std::io::Cursor;

#[test]
fn query_intervals_follow_alignment_direction() {
    for offset in 0..=22 {
        for length in 0..=22 - offset {
            assert_eq!(
                query_interval(Strand::Positive, 88, 110, offset, length),
                (88 + offset, 88 + offset + length)
            );
            assert_eq!(
                query_interval(Strand::Negative, 88, 110, offset, length),
                (110 - offset - length, 110 - offset)
            );
        }
    }
}

fn call_fixture(chunk_size: usize) -> Vec<Vec<String>> {
    let fixture = include_str!("../../test/reverse_strand_indels.maf");
    let mut reader = MAFReader::new(Cursor::new(fixture.as_bytes())).unwrap();
    let mut output = Vec::new();
    call_var_maf(
        &mut reader,
        None,
        &mut output,
        true,
        false,
        0,
        false,
        None,
        None,
        None,
        Some(chunk_size),
    )
    .unwrap();
    let text = String::from_utf8(output).unwrap();
    let mut records: Vec<Vec<String>> = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').map(str::to_owned).collect())
        .collect();
    records.sort();
    records
}

#[test]
fn maf_calls_preserve_alleles_and_report_forward_query_coordinates() {
    let records = call_fixture(1_000_000);
    let expected = [
        ("chr1", "61", "C", "CTTGACCGATGCTAAC", "1|1:chr1A@50@65@P"),
        ("chr1", "66", "GTACGA", "G", "1|1:chr1A@70@70@P"),
        ("chr10", "15", "A", "ATTGG", "1|1:chr10A@71@75@N"),
        ("chr4", "36", "C", "CTTGGA", "1|1:chr4A@99@104@N"),
        ("chr4", "42", "TACG", "T", "1|1:chr4A@93@93@N"),
        ("chrS", "13", "G", "A", "1|1:chrSA@77@N"),
    ];
    assert_eq!(records.len(), expected.len());
    for (chrom, pos, reference, alternate, sample) in expected {
        let record = records
            .iter()
            .find(|row| row[0] == chrom && row[1] == pos)
            .unwrap();
        assert_eq!(record[3], reference);
        assert_eq!(record[4], alternate);
        assert_eq!(record[9], sample);
    }
}

#[test]
fn maf_query_coordinates_are_independent_of_chunk_size() {
    assert_eq!(call_fixture(7), call_fixture(1_000_000));
}
