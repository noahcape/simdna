use std::arch::aarch64::uint8x16_t;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use simdna::seed::*;

pub fn file_scan_simd(c: &mut Criterion) {
    use std::{
        fs::File,
        io::{BufRead, BufReader},
    };

    let pattern: uint8x16_t = unsafe { SIMDna::load_pattern(b"CAGAGC") };

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

                while start < seq.len() - 16 {
                    let seq: uint8x16_t =
                        unsafe { SIMDna::load_ref(seq[start..start + 16].as_bytes()) };
                    unsafe {
                        pattern
                            .shuffle_bytes(seq)
                            .shift_lanes()
                            .count_ones()
                            .locate_seeds(3)
                            .get_seeds(start)
                    };

                    start += 10;
                }

                if start < seq.len() {
                    let seq: uint8x16_t =
                        unsafe { SIMDna::load_ref(seq[seq.len() - 16..seq.len()].as_bytes()) };
                    unsafe {
                        pattern
                            .shuffle_bytes(seq)
                            .shift_lanes()
                            .count_ones()
                            .locate_seeds(3)
                            .get_seeds(start)
                    };
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

// in order to replicate the current scanning with ANTISEQUENCE
fn hamming(a: &[u8], b: &[u8], threshold: usize) -> Option<usize> {
    if a.len() != b.len() {
        return None;
    }

    let a_ptr = a.as_ptr();
    let b_ptr = b.as_ptr();
    let n = a.len();
    let mut res = 0;
    let mut i = 0;

    unsafe {
        while i < (n / 8) * 8 {
            let a_word = std::ptr::read_unaligned(a_ptr.add(i) as *const u64);
            let b_word = std::ptr::read_unaligned(b_ptr.add(i) as *const u64);

            let xor = a_word ^ b_word;
            let or1 = xor | (xor >> 1);
            let or2 = or1 | (or1 >> 2);
            let or3 = or2 | (or2 >> 4);
            let mask = or3 & 0x0101010101010101u64;
            res += mask.count_ones() as usize;

            i += 8;
        }

        while i < n {
            res += (*a_ptr.add(i) != *b_ptr.add(i)) as usize;
            i += 1;
        }
    }

    let matches = n - res;

    if matches >= threshold {
        Some(matches)
    } else {
        None
    }
}
