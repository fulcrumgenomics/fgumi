//! `OrderedBytesSingle<T>` requires `T: HeapSize + Ordered`. A type that
//! impls `Ordered` but NOT `HeapSize` must fail to compile.

use fgumi_lib::pipeline::core::item::Ordered;
use fgumi_lib::pipeline::core::outputs::OrderedBytesSingle;

struct OrderedNoHeap {
    serial: u64,
}

impl Ordered for OrderedNoHeap {
    fn ordinal(&self) -> u64 {
        self.serial
    }
}

fn _instantiate() -> std::marker::PhantomData<OrderedBytesSingle<OrderedNoHeap>> {
    // Bound `T: HeapSize` is unsatisfied for `OrderedNoHeap`.
    std::marker::PhantomData
}

fn main() {}
