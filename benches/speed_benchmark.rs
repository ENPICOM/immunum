use criterion::{criterion_group, criterion_main, BatchSize, Criterion, Throughput};
use immunum::types::Chain;
use immunum::{read_input, Annotator, Scheme};
use std::time::Duration;

const FASTA: &str = "fixtures/SRR_HKL_10k.fasta";

fn bench_annotate(c: &mut Criterion) {
    let records = read_input(Some(FASTA)).expect("Failed to read FASTA");
    let sequences: Vec<String> = records.into_iter().map(|r| r.sequence).collect();

    let chains = [Chain::IGH, Chain::IGK, Chain::IGL];
    let annotator = Annotator::new(&chains, Scheme::IMGT).expect("Failed to create annotator");

    let mut group = c.benchmark_group("annotate");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(30));
    group.throughput(Throughput::Elements(sequences.len() as u64));

    group.bench_function("SRR_HKL_10k", |b| {
        b.iter_batched(
            || sequences.clone(),
            |seqs| {
                for seq in &seqs {
                    let _ = annotator.number(seq);
                }
            },
            BatchSize::LargeInput,
        )
    });

    group.finish();
}

criterion_group!(benches, bench_annotate);
criterion_main!(benches);
