use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_pgen::Pgen;
use rsomics_plink_ibc::ibc;
use std::hint::black_box;
use std::path::PathBuf;

fn bench_ibc(c: &mut Criterion) {
    // Set RSOMICS_IBC_BFILE to a representative-large fileset prefix; otherwise
    // the in-repo golden is used.
    let prefix = std::env::var("RSOMICS_IBC_BFILE")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/small"));
    let pgen = Pgen::load(&prefix).expect("load fileset");

    c.bench_function("ibc", |b| {
        b.iter(|| black_box(ibc(black_box(&pgen))));
    });
}

criterion_group!(benches, bench_ibc);
criterion_main!(benches);
