use std::arch::aarch64::uint8x16_t;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use simdna::{hamming, simdna::*, fallback::Patterns};

pub fn file_scan_simd(c: &mut Criterion) {
    use std::{
        fs::File,
        io::{BufRead, BufReader},
    };

    let pattern_seq = b"CAGAGC";

    let pattern: uint8x16_t = unsafe { SIMDna::load_pattern(pattern_seq) };

    let file_path = "./seq_lib.txt";

    let file = File::open(file_path).expect("no such file");
    let buf = BufReader::new(file);
    let lines: Vec<String> = buf
        .lines()
        .map(|l| l.expect("Could not parse line"))
        .collect();

    c.bench_function("scan file with SIMD", |b| {
        b.iter(|| {
            for seq in &lines {
                let mut start = 0;

                if let Some(idx) = unsafe { pattern.locate(seq.as_bytes(), 3, &mut start) } {
                    hamming(pattern_seq, seq[idx..idx + pattern_seq.len()].as_bytes(), 5);
                }
            }
        })
    });
}

pub fn file_scan(c: &mut Criterion) {
    use std::{
        fs::File,
        io::{BufRead, BufReader},
    };

    let pattern = Patterns::new(b"CAGAGC", 13);

    let file_path = "./seq_lib.txt";

    let file = File::open(file_path).expect("no such file");
    let buf = BufReader::new(file);
    let lines: Vec<String> = buf
        .lines()
        .map(|l| l.expect("Could not parse line"))
        .collect();

    c.bench_function("scan file", |b| {
        b.iter(|| {
            for seq in &lines {
                pattern.seed(black_box(seq.as_bytes()));
            }
        })
    });
}

pub fn file_scan_edit_dist(c: &mut Criterion) {
    use std::{
        fs::File,
        io::{BufRead, BufReader},
    };

    let pattern = b"CAGAGC";

    let file_path = "./seq_lib.txt";

    let file = File::open(file_path).expect("no such file");
    let buf = BufReader::new(file);
    let lines: Vec<String> = buf
        .lines()
        .map(|l| l.expect("Could not parse line"))
        .collect();

    c.bench_function("scan file edit dist", |b| {
        b.iter(|| {
            let mut best_match = None;

            for seq in &lines {
                for (i, overlap) in seq.as_bytes()[..].windows(pattern.len()).enumerate() {
                    if let Some(matches) = hamming(overlap, pattern, 5) {
                        if let Some((best_matches, _, _)) = best_match {
                            if matches <= best_matches {
                                continue;
                            }
                        }

                        best_match = Some((matches, i, i + pattern.len()));
                    }
                }
            }
        })
    });
}

criterion_group!(benches, file_scan_simd, file_scan, file_scan_edit_dist);

criterion_main!(benches);
