pub mod seed;

use std::arch::aarch64::uint8x16_t;

use crate::seed::SIMDna;

#[inline]
pub fn main() {
    let pattern_seq = b"CAGAGC";
    let pattern: uint8x16_t = unsafe { SIMDna::load_pattern(pattern_seq) };

    let seq = "TCAGATTCTCCCCGGATTTAATTGAATTTT";

    let mut start = 0;

    loop {
        if let Some(idx) = unsafe { pattern.locate(seq.as_bytes(), 3, &mut start) } {
            if let Some(_) = hamming(pattern_seq, seq[idx..idx + pattern_seq.len()].as_bytes(), 5) {
                println!("Match at: {}", idx);
                break;
            }

            start = idx + 1;
        } else {
            break;
        }
    }
}

// Hamming distance algorithm @daniel_c0db0t in ANTISEQUENCE library
pub fn hamming(a: &[u8], b: &[u8], threshold: usize) -> Option<usize> {
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
