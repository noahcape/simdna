pub trait SIMDna {
    fn block_size() -> usize;

    /// used scheme c[i] = a[b[i]] where c is returned `Self`, b are `indices` and `a` is being indexed
    unsafe fn shuffle_bytes(self, indices: Self) -> Self;

    /// sets all lanes whose byte population is greater than threshold to all ones
    unsafe fn fill_seed_lanes(self, threshold: u8) -> Self;

    /// locate the first instance of a seed
    /// `offset` indicates lane[0]'s distance from beginning of sequence
    /// `len` is length of input seq
    unsafe fn find(self, offset: usize, len: usize) -> Option<usize>;

    /// shifting lanes 8 times to align down sequence kmer instances with start instance
    unsafe fn shift_lanes(self) -> Self;

    /// Creates a new compressed nucleotide kmer SIMD vector, `slice` should be at least a 17 base dna sequence
    unsafe fn load_ref(slice: &[u8]) -> Self;

    /// Creates a new finger print pattern representation
    unsafe fn load_pattern(slice: &[u8]) -> Self;

    /// Turn array of u8 into SIMD vector
    unsafe fn load(slice: &[u8; 16]) -> Self;

    /// Run methods to find first acceptable seed index, threshold is passed to `fill_seed_lanes`
    /// This will only scan ref_ up until first seed is found, if it is not acceptable run again with `start` as the return index + 1
    unsafe fn locate(self, ref_: &[u8], threshold: u8, start: &mut usize) -> Option<usize>;
}
