#[derive(Clone, Debug)]
pub struct Patterns {
    patterns: [u8; 16],
    threshold: u32,
}

impl Patterns {
    pub fn new(a: &[u8]) -> Self {
        let mask = 0x06;
        let mut patterns = [0u8; 16];

        let n = if a.len() > 8 { 8 } else { a.len() };

        for (i, fp) in a[0..=n].windows(2).enumerate() {
            let idx = (fp[0] & mask) << 1 | (fp[1] & mask) >> 1;
            patterns[idx as usize] |= 1 << i
        }

        Self {
            patterns,
            threshold: (a.len() / 2) as u32,
        }
    }

    pub fn seed(&self, ref_: &[u8], _identity: f64) -> Vec<usize> {
        let mut seeds: Vec<usize> = vec![];

        if ref_.len() <= 16 {
            let c = self.pshufb(ref_);

            seeds.append(&mut self.extract_seeds(&c, 0))
        } else {
            let mut start = 0;

            while start < ref_.len() - 16 {
                let c = self.pshufb(&ref_[start..=start + 16]);

                seeds.append(&mut self.extract_seeds(&c, start));

                start += 16;
            }
        }

        seeds
    }

    fn pshufb(&self, ref_: &[u8]) -> [u8; 16] {
        let mut c: [u8; 16] = [0u8; 16];

        for (i, kmer) in ref_[..ref_.len() - 1].windows(2).enumerate() {
            let idx = ((kmer[0] & 0x6) << 1) | ((kmer[1] & 0x6) >> 1);
            let fp = self.patterns[idx as usize] as u128;

            let to = if 8 > i { i } else { 8 };

            for j in 0..to {
                let mask = 1 << j;
                c[i - j] |= (fp & mask) as u8;
            }
        }

        c
    }

    fn extract_seeds(&self, c: &[u8; 16], offset: usize) -> Vec<usize> {
        let mut pos: Vec<usize> = vec![];

        for i in 0..c.len() {
            let idx = i + offset;
            if c[i].count_ones() >= self.threshold {
                pos.push(idx)
            }
        }

        pos
    }
}
