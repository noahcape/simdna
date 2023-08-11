use std::arch::aarch64::uint8x16_t;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use simdna::{hamming, seed::*};

pub fn simd_byte_pattern_scan(c: &mut Criterion) {
    let pattern_seq = b"CAGAGC";

    let pattern: PackedByte<uint8x16_t> = unsafe { SIMDByte::load_pattern(pattern_seq) };

    c.bench_function("pattern scan with simd byte", |b| {
        b.iter(|| {
            let seq = "TCAGATTCTCCCCGGATTTAATTGAATTTT";

            let mut start = 0;

            loop {
                if let Some(idx) = unsafe {
                    pattern.locate(
                        black_box(seq.as_bytes()),
                        black_box(5),
                        black_box(&mut start),
                    )
                } {
                    if idx + pattern_seq.len() <= seq.len() {
                        if let Some(_) = hamming(
                            black_box(pattern_seq),
                            black_box(seq[idx..idx + pattern_seq.len()].as_bytes()),
                            black_box(5),
                        ) {
                            break;
                        }
                    }

                    start = idx + 1;
                } else {
                    break;
                }
            }
        })
    });
}

pub fn simd_pattern_scan(c: &mut Criterion) {
    let pattern_seq = b"CAGAGC";

    let pattern: uint8x16_t = unsafe { SIMDna::load_pattern(pattern_seq) };

    c.bench_function("pattern scan with simd", |b| {
        b.iter(|| {
            let seq = "TCAGATTCTCCCCGGATTTAATTGAATTTT";

            let mut start = 0;

            loop {
                if let Some(idx) = unsafe {
                    pattern.locate(
                        black_box(seq.as_bytes()),
                        black_box(3),
                        black_box(&mut start),
                    )
                } {
                    if idx + pattern_seq.len() <= seq.len() {
                        if let Some(_) = hamming(
                            black_box(pattern_seq),
                            black_box(seq[idx..idx + pattern_seq.len()].as_bytes()),
                            black_box(5),
                        ) {
                            break;
                        }
                    }

                    start = idx + 1;
                } else {
                    break;
                }
            }
        })
    });
}

pub fn hamming_patter_scan(c: &mut Criterion) {
    let pattern = b"CAGAGC";

    c.bench_function("pattern scan with hamming dist", |b| {
        b.iter(|| {
            let mut best_match = None;

            let seq = "TCAGATTCTCCCCGGATTTAATTGAATTTT";

            for (i, overlap) in seq.as_bytes()[..].windows(pattern.len()).enumerate() {
                if let Some(matches) = hamming(black_box(overlap), black_box(pattern), black_box(5))
                {
                    if let Some((best_matches, _, _)) = best_match {
                        if matches <= best_matches {
                            continue;
                        }
                    }

                    best_match = Some((matches, i, i + pattern.len()));
                }
            }
        })
    });
}

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

                loop {
                    if let Some(idx) = unsafe { pattern.locate(seq.as_bytes(), 3, &mut start) } {
                        if idx + pattern_seq.len() <= seq.len() {
                            if let Some(_) = hamming(
                                black_box(pattern_seq),
                                black_box(seq[idx..idx + pattern_seq.len()].as_bytes()),
                                black_box(5),
                            ) {
                                break;
                            }
                        }

                        start = idx + 1;
                    } else {
                        break;
                    }
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

criterion_group!(
    benches,
    simd_byte_pattern_scan,
    simd_pattern_scan,
    hamming_patter_scan,
    file_scan_simd,
    file_scan,
    file_scan_edit_dist
);

criterion_main!(benches);
