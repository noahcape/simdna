// First attempt at seeding without SIMD
#[derive(Clone, Debug)]
pub struct Patterns {
    patterns: [u8; 16],
    threshold: u32,
    len: usize,
}

impl Patterns {
    pub fn new(a: &[u8], threshold: usize) -> Self {
        let mask = 0x06;
        let mut patterns = [0u8; 16];

        let n = if a.len() > 8 { 8 } else { a.len() };

        for (i, fp) in a[0..n].windows(2).enumerate() {
            let idx = (fp[0] & mask) << 1 | (fp[1] & mask) >> 1;
            patterns[idx as usize] |= 1 << i
        }

        Self {
            patterns,
            threshold: (n * a.len() / threshold / 2) as u32,
            len: a.len(),
        }
    }

    pub fn seed(&self, ref_: &[u8]) -> Vec<usize> {
        if ref_.len() <= 16 {
            return self.extract_seeds(ref_, 0, ref_.len());
        }

        let mut start = 0;
        let mut seeds: Vec<usize> = vec![];

        while start < ref_.len() - 16 {
            let mut c = self.extract_seeds(&ref_[start..=start + 16], start, ref_.len());

            seeds.append(&mut c);

            start += 16 - self.len;
        }

        if start < ref_.len() {
            let mut c = self.extract_seeds(&ref_[start..], start, ref_.len());

            seeds.append(&mut c);
        }

        seeds
    }

    fn extract_seeds(&self, ref_: &[u8], offset: usize, len: usize) -> Vec<usize> {
        let mut c: [u8; 16] = [0u8; 16];
        let mut pos: Vec<usize> = vec![];

        for (i, kmer) in ref_[..ref_.len()].windows(2).enumerate() {
            let idx = ((kmer[0] & 0x6) << 1) | ((kmer[1] & 0x6) >> 1);
            let fp = self.patterns[idx as usize];

            let to = if 8 > i { i + 1 } else { 8 };

            for j in 0..to {
                let mask = 1 << j;
                if fp & mask == mask {
                    c[i - j] |= fp & mask;
                }
            }

            if i >= to && c[i - to].count_ones() >= self.threshold && i - to + self.len <= len {
                pos.push(i - to + offset)
            }
        }

        for (i, elm) in c
            .iter()
            .enumerate()
            .take(ref_.len() - 1)
            .skip(ref_.len() - 9)
        {
            if elm.count_ones() >= self.threshold && i + self.len <= len {
                pos.push(i + offset)
            }
        }

        pos
    }
}

#[test]
fn test_seed() {
    let fp = Patterns::new(b"CAGAGC", 6);

    let seeds = fp.seed(b"TATAAGGCCTGTCTCTTATACACATCTCCGAGCCCA");

    assert_eq!(vec![27], seeds);
}