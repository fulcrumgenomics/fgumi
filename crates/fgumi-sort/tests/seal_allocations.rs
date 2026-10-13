//! A parallel template seal allocates the chunk's output tuples and the
//! radix's own scratch, and nothing more: the merge writes the tuples
//! directly, with no intermediate sorted-refs buffer and no copy back; and it
//! frees the refs, so the sealed run holds one run-sized buffer. Its own
//! binary, because it installs a counting global allocator (`dhat`, in testing
//! mode).
#![cfg(not(loom))]

use fgumi_sort::{
    TemplateKey24, TemplateRecordRef, radix_sort_template_refs, radix_sorts_in_parallel,
    sort_template_refs_into_records,
};
use rayon::prelude::*;

#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

#[test]
fn seal_allocates_only_the_output_tuples() {
    const THREADS: usize = 4;
    let n = 512 * 1024;
    assert!(radix_sorts_in_parallel(n, THREADS), "the fixture takes the parallel path");
    let refs: Vec<TemplateRecordRef<TemplateKey24>> = (0..n)
        .map(|i| {
            let x = (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
            TemplateRecordRef {
                key: TemplateKey24 { primary: x >> 40, secondary: x & 0xFF, name_hash_upper: x },
                offset: i as u64 * 16,
                len: 100,
                padding: 0,
            }
        })
        .collect();
    let pool = rayon::ThreadPoolBuilder::new().num_threads(THREADS).build().unwrap();
    let _profiler = dhat::Profiler::builder().testing().build();

    // The radix's own scratch, sorting the same per-thread chunks the seal does.
    let mut chunks = refs.clone();
    let before = dhat::HeapStats::get().total_bytes;
    pool.install(|| {
        chunks.par_chunks_mut(n.div_ceil(THREADS)).for_each(radix_sort_template_refs);
    });
    let scratch = dhat::HeapStats::get().total_bytes - before;
    drop(chunks);

    // Re-allocate the refs under the profiler, so their free is counted.
    let refs = refs.clone();
    let refs_bytes = (refs.capacity() * size_of::<TemplateRecordRef<TemplateKey24>>()) as u64;
    let live_before = dhat::HeapStats::get().curr_bytes as u64;
    let before = dhat::HeapStats::get().total_bytes;
    let records = pool.install(|| sort_template_refs_into_records(refs, THREADS));
    let allocated = dhat::HeapStats::get().total_bytes - before;
    let live_after = dhat::HeapStats::get().curr_bytes as u64;

    let output = (records.capacity() * size_of::<(TemplateKey24, u64, u32)>()) as u64;
    assert_eq!(records.len(), n);
    assert!(allocated >= output, "the output tuples are allocated ({allocated} < {output})");
    let slack = 64 * 1024;
    assert!(
        allocated <= output + scratch + slack,
        "the seal allocated {allocated} bytes: output {output} + radix scratch {scratch} + \
         {slack} slack — an intermediate sorted-refs buffer would add {}",
        n * size_of::<TemplateRecordRef<TemplateKey24>>()
    );
    // The pool's worker threads may keep a few small allocations alive (thread
    // locals, job bookkeeping), hence the slack; a retained ref buffer would
    // leave all `refs_bytes` live.
    let expected = live_before - refs_bytes + output;
    assert!(
        live_after.abs_diff(expected) <= slack,
        "after the seal only the output tuples are live (got {live_after}, expected \
         {expected} ± {slack}): the refs ({refs_bytes} bytes) are freed"
    );
}
