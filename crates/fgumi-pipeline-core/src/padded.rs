//! Cache-line padding shared by the crate's hot counters.

/// A value alone on its own 128-byte line: 128 rather than 64 because Apple
/// Silicon and some `x86_64` prefetchers pair adjacent 64-byte lines, so a
/// 64-byte stride can still share a coherence unit. The one padding type of
/// the crate: the event-count, the admission caps, the liveness counters and
/// the wake slots use it, so their hot atomics never false-share.
#[repr(align(128))]
#[derive(Debug, Default)]
pub(crate) struct Padded<T>(pub(crate) T);

const _: () = assert!(std::mem::size_of::<Padded<std::sync::atomic::AtomicUsize>>() == 128);
const _: () = assert!(std::mem::size_of::<Padded<std::sync::atomic::AtomicU64>>() == 128);
