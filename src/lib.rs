pub mod simdna;

// SIMD instructions
#[cfg(feature = "simd_neon")]
pub mod neon;

#[cfg(feature = "simd_avx2")]
pub mod avx2;

pub mod fallback;
pub mod utils;

use utils::hamming;
use std::arch::aarch64::uint8x16_t;

use crate::simdna::SIMDna;

#[inline]
pub fn main() {
    let pattern_seq = b"CAGAGC";
    let pattern: uint8x16_t = unsafe { SIMDna::load_pattern(pattern_seq) };

    let seq = "TCAGATTCTCCCCGGATTTAATCAGAGCTGAATTTT";

    let mut start = 0;

    loop {
        if let Some(idx) = unsafe { pattern.locate(seq.as_bytes(), 3, &mut start) } {
            if let Some(_) = hamming(pattern_seq, seq[idx..idx + pattern_seq.len()].as_bytes(), 5) {
                println!("Match at: {}", idx);
                break;
            }

            start = idx + 1;
        }
    }
}
