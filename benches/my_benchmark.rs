use criterion::{black_box, criterion_group, criterion_main, Criterion};
use simdna::{
    seed::Patterns
};

pub fn file_scan(c: &mut Criterion) {
    use std::{
        fs::File,
        io::{BufRead, BufReader},
    };

    let pattern = Patterns::new(b"AGCTCATAAGTGAG", 13);

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
criterion_group!(
    benches,
    file_scan,
);

criterion_main!(benches);
