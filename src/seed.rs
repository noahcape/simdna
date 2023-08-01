#![allow(unused_imports)]

#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
#[cfg(target_arch = "aarch64")]
use std::{arch::aarch64::*, mem::transmute};

static SHIFT_LANE_LEFT: [[u8; 16]; 8] = [
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0],
    [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 0],
    [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 0, 0],
    [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 0, 0, 0],
    [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 0, 0, 0, 0],
    [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 0, 0, 0, 0, 0],
    [7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 0, 0, 0, 0, 0, 0],
];

pub trait SIMDna {
    /// used scheme c[i] = a[b[i]] where c is returned `Self`, b are `indices` and `a` is being indexed
    fn shuffle_bytes(self, indices: Self) -> Self;

    /// count set bits in each lane (byte)
    fn count_ones(self) -> Self;

    fn locate_seeds(self, threshold: u8) -> Self;
    fn get_seeds(self) -> u8;
    unsafe fn shift_lanes(self) -> Self;

    /// Creates a new compressed nucleotide kmer SIMD vector, `slice` should be 17 base dna sequence
    unsafe fn load_ref(slice: &[u8]) -> Self;

    /// Creates a new finger print pattern representation
    unsafe fn load_pattern(slice: &[u8]) -> Self;
}

impl SIMDna for uint8x16_t {
    #[inline]
    fn shuffle_bytes(self, indices: Self) -> Self {
        unsafe { vqtbl1q_u8(self, indices) }
    }

    #[inline]
    fn count_ones(self) -> Self {
        unsafe { vcntq_u8(self) }
    }

    #[inline]
    fn locate_seeds(self, threshold: u8) -> Self {
        unsafe { vcgeq_u8(self, transmute([threshold; 16])) }
    }

    #[inline]
    fn get_seeds(self) -> u8 {
        // this would be a way to traverse lanes
        unsafe { vgetq_lane_u8::<5>(self) }
    }

    #[inline]
    unsafe fn shift_lanes(self) -> Self {
        let mut res: Self = transmute([0u8; 16]);

        for i in 0..8 {
            res = vorrq_u8(
                res,
                vandq_u8(self, transmute([(1 << i) as u8; 16]))
                    .shuffle_bytes(transmute(SHIFT_LANE_LEFT[i])),
            );
        }

        res
    }

    #[inline]
    unsafe fn load_ref(slice: &[u8]) -> uint8x16_t {
        let mut b = [0u8; 16];

        for (i, kmer) in slice.windows(2).enumerate() {
            b[i] = ((kmer[0] & 0x6) << 1) | ((kmer[1] & 0x6) >> 1);
        }

        transmute(b)
    }

    #[inline]
    unsafe fn load_pattern(slice: &[u8]) -> Self {
        let mask = 0x06;
        let mut patterns = [0u8; 16];

        let n = if slice.len() > 8 { 8 } else { slice.len() };

        for (i, fp) in slice[0..n].windows(2).enumerate() {
            let idx = (fp[0] & mask) << 1 | (fp[1] & mask) >> 1;
            patterns[idx as usize] |= 1 << i
        }

        transmute(patterns)
    }
}

#[test]
fn simd_instr() {
    unsafe {
        let pattern_vec: uint8x16_t = SIMDna::load_pattern(b"CAGAGC");
        let ref_vec: uint8x16_t = SIMDna::load_ref(b"TTTTTCAGAGCTTTTTT");
        let c = pattern_vec.shuffle_bytes(ref_vec);

        println!("{:?}", c.shift_lanes());
        println!("{:?}", c.shift_lanes().count_ones());
        println!("{:?}", c.shift_lanes().count_ones().locate_seeds(3));
        println!(
            "{:?}",
            c.shift_lanes().count_ones().locate_seeds(3).get_seeds()
        );
    }
}

/// First attempt at seeding without SIMD
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
