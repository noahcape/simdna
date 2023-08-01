use criterion::{black_box, criterion_group, criterion_main, Criterion};
use simdna::seed::Patterns;

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
            for seq in &lines {
                for overlap in seq.as_bytes()[..].windows(pattern.len()) {
                    hamming(black_box(overlap), black_box(pattern), black_box(13));
                }
            }
        })
    });
}

criterion_group!(benches, file_scan, file_scan_edit_dist);

criterion_main!(benches);

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
