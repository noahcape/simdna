use std::arch::aarch64::*;

use std::{mem::size_of, ptr};

use crate::simdna::SIMDna;

type SIMD = uint8x16_t;

static SHIFT_LANE_LEFT: [&[u8; 16]; 8] = [
    &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
    &[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0],
    &[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 0],
    &[3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 0, 0],
    &[4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 0, 0, 0],
    &[5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 0, 0, 0, 0],
    &[6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 0, 0, 0, 0, 0],
    &[7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 0, 0, 0, 0, 0, 0],
];

static LANE_AND: [&[u8; 16]; 8] = [
    &[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    &[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0],
    &[4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0],
    &[8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 0, 0, 0],
    &[16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 0, 0, 0, 0],
    &[32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 0, 0, 0, 0, 0],
    &[64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 0, 0, 0, 0, 0, 0],
    &[
        128, 128, 128, 128, 128, 128, 128, 128, 128, 0, 0, 0, 0, 0, 0, 0,
    ],
];

impl SIMDna for uint8x16_t {
    #[inline]
    fn block_size() -> usize {
        16
    }

    #[target_feature(enable = "neon")]
    #[inline]
    unsafe fn shuffle_bytes(self, indices: Self) -> Self {
        vqtbl1q_u8(self, indices)
    }

    #[target_feature(enable = "neon")]
    #[inline]
    unsafe fn fill_seed_lanes(self, threshold: u8) -> Self {
        vcgeq_u8(vcntq_u8(self), SIMDna::load(&[threshold; 16]))
    }

    #[target_feature(enable = "neon")]
    #[inline]
    unsafe fn find(self, offset: usize, len: usize) -> Option<usize> {
        if vmaxvq_u8(self) != 255 {
            return None;
        }

        let arr = ptr::addr_of!(self) as *const u8;

        let mut i = 0;
        while i < len {
            if *arr.add(i * size_of::<u8>()) == 255 {
                return Some(i + offset);
            }

            i += 1;
        }

        None
    }

    #[target_feature(enable = "neon")]
    #[inline]
    unsafe fn shift_lanes(self) -> Self {
        let mut res: Self = SIMDna::load(&[0u8; 16]);

        for i in 0..8 {
            res = vorrq_u8(
                res,
                vandq_u8(self, SIMDna::load(LANE_AND[i]))
                    .shuffle_bytes(SIMDna::load(SHIFT_LANE_LEFT[i])),
            );
        }

        res
    }

    #[target_feature(enable = "neon")]
    #[inline]
    unsafe fn locate(self, ref_: &[u8], threshold: u8, start: &mut usize) -> Option<usize> {
        let block_size = SIMD::block_size();

        while *start < ref_.len() - block_size {
            let seq: SIMD = SIMDna::load_ref(&ref_[*start..=*start + block_size]);
            let idx = self
                .shuffle_bytes(seq)
                .shift_lanes()
                .fill_seed_lanes(threshold)
                .find(*start, block_size);

            if idx.is_some() {
                return idx;
            }

            *start += 10;
        }

        if *start >= ref_.len() - block_size {
            let seq_vec: SIMD = SIMDna::load_ref(&ref_[*start..]);
            let idx = self
                .shuffle_bytes(seq_vec)
                .shift_lanes()
                .fill_seed_lanes(3)
                .find(*start, ref_.len() - *start);

            if idx.is_some() {
                return idx;
            }
        }

        None
    }

    #[target_feature(enable = "neon")]
    #[inline]
    unsafe fn load_ref(slice: &[u8]) -> SIMD {
        let mut b = [0u8; 16];

        for (i, kmer) in slice.windows(2).enumerate() {
            b[i] = ((kmer[0] & 0x6) << 1) | ((kmer[1] & 0x6) >> 1);
        }

        SIMDna::load(&b)
    }

    #[target_feature(enable = "neon")]
    #[inline]
    unsafe fn load_pattern(slice: &[u8]) -> Self {
        let mask = 0x06;
        let mut patterns = [0u8; 16];

        let n = if slice.len() > 8 { 8 } else { slice.len() };

        for (i, fp) in slice[0..n].windows(2).enumerate() {
            let idx = (fp[0] & mask) << 1 | (fp[1] & mask) >> 1;
            patterns[idx as usize] |= 1 << i
        }

        SIMDna::load(&patterns)
    }

    #[target_feature(enable = "neon")]
    #[inline]
    unsafe fn load(slice: &[u8; 16]) -> Self {
        vld1q_u8(slice.as_ptr() as *const u8)
    }
}


// TODO better tests
#[test]
fn simd_instr() {
    unsafe {
        let pattern_vec: uint8x16_t = SIMDna::load_pattern(b"CAGAGC");
        let ref_vec: uint8x16_t = SIMDna::load_ref(b"TCAGAGCTTTTTTTTTT");
        let c = pattern_vec.shuffle_bytes(ref_vec);

        println!("{:?}", c.shift_lanes().fill_seed_lanes(3).find(0, 16));
    }
}
