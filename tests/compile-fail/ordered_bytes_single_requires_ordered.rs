//! `OrderedBytesSingle<T>` requires `T: HeapSize + Ordered`. A type that
//! impls `HeapSize` but NOT `Ordered` must fail to compile.

use fgumi_lib::pipeline::core::item::HeapSize;
use fgumi_lib::pipeline::core::outputs::OrderedBytesSingle;

struct HeapNoOrdered {
    bytes: Vec<u8>,
}

impl HeapSize for HeapNoOrdered {
    fn heap_size(&self) -> usize {
        self.bytes.len()
    }
}

fn _instantiate() -> std::marker::PhantomData<OrderedBytesSingle<HeapNoOrdered>> {
    // Bound `T: Ordered` is unsatisfied for `HeapNoOrdered`.
    std::marker::PhantomData
}

fn main() {}
