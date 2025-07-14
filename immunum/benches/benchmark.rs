use criterion::{criterion_group, criterion_main, Criterion};
use immunum::numbering::number_sequence;
use immunum::types::{Chain, Scheme};
use std::hint::black_box;

fn benchmark_numbering_short_sequence(c: &mut Criterion) {
    let sequence = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGW";
    let scheme = Scheme::IMGT;
    let chains = vec![Chain::IGH];

    c.bench_function("number_sequence_short", |b| {
        b.iter(|| number_sequence(black_box(sequence), black_box(&scheme), black_box(&chains)))
    });
}

fn benchmark_numbering_long_sequence(c: &mut Criterion) {
    let sequence = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCARGGYDILTGYYFDYWGQGTLVTVSS";
    let scheme = Scheme::IMGT;
    let chains = vec![Chain::IGH];

    c.bench_function("number_sequence_long", |b| {
        b.iter(|| number_sequence(black_box(sequence), black_box(&scheme), black_box(&chains)))
    });
}

fn benchmark_numbering_different_schemes(c: &mut Criterion) {
    let sequence = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGW";
    let chains = vec![Chain::IGH];

    c.bench_function("number_sequence_imgt", |b| {
        b.iter(|| {
            number_sequence(
                black_box(sequence),
                black_box(&Scheme::IMGT),
                black_box(&chains),
            )
        })
    });

    c.bench_function("number_sequence_kabat", |b| {
        b.iter(|| {
            number_sequence(
                black_box(sequence),
                black_box(&Scheme::KABAT),
                black_box(&chains),
            )
        })
    });
}

fn benchmark_numbering_different_chains(c: &mut Criterion) {
    let sequence = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGW";
    let scheme = Scheme::IMGT;

    c.bench_function("number_sequence_heavy_chain", |b| {
        b.iter(|| {
            number_sequence(
                black_box(sequence),
                black_box(&scheme),
                black_box(&vec![Chain::IGH]),
            )
        })
    });

    c.bench_function("number_sequence_kappa_chain", |b| {
        b.iter(|| {
            number_sequence(
                black_box(sequence),
                black_box(&scheme),
                black_box(&vec![Chain::IGK]),
            )
        })
    });

    c.bench_function("number_sequence_lambda_chain", |b| {
        b.iter(|| {
            number_sequence(
                black_box(sequence),
                black_box(&scheme),
                black_box(&vec![Chain::IGL]),
            )
        })
    });
}

criterion_group!(
    benches,
    benchmark_numbering_short_sequence,
    benchmark_numbering_long_sequence,
    benchmark_numbering_different_schemes,
    benchmark_numbering_different_chains
);
criterion_main!(benches);
