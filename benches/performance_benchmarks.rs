use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use immunum::{
    annotator::Annotator,
    types::{Chain, Scheme},
};
use std::{hint::black_box, time::Duration};

/// Generate sequences of specific lengths for testing algorithmic complexity
fn create_sequence_of_length(length: usize) -> String {
    // Base heavy chain sequence
    let base = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARWGGRGSYAMDYWGQGTLVTVSS";

    if length <= base.len() {
        base[..length].to_string()
    } else {
        // For longer sequences, repeat and add realistic amino acids
        let mut sequence = base.to_string();
        let extension = "AKTTAPSVYPLAPLSPSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKS";

        while sequence.len() < length {
            let remaining = length - sequence.len();
            if remaining >= extension.len() {
                sequence.push_str(extension);
            } else {
                sequence.push_str(&extension[..remaining]);
            }
        }
        sequence
    }
}

/// Create test sequences with realistic length distribution
fn create_realistic_test_sequences(count: usize) -> Vec<String> {
    // Realistic antibody sequence lengths based on biological data
    let realistic_sequences = vec![
        // Short variable regions (VH/VL domains only) - 100-150 aa
        create_sequence_of_length(110),
        create_sequence_of_length(105),
        create_sequence_of_length(130),
        create_sequence_of_length(125),
        // Medium sequences (V + partial constant) - 200-300 aa
        create_sequence_of_length(220),
        create_sequence_of_length(250),
        create_sequence_of_length(280),
        // Long sequences (full heavy chain) - 400-500+ aa
        create_sequence_of_length(450),
        create_sequence_of_length(480),
        create_sequence_of_length(520),
    ];

    (0..count)
        .map(|i| realistic_sequences[i % realistic_sequences.len()].clone())
        .collect()
}

/// Benchmark single sequence processing
fn benchmark_single_sequence(c: &mut Criterion) {
    let annotator = Annotator::new(
        Scheme::IMGT,
        vec![Chain::IGH, Chain::IGK, Chain::IGL],
        None,
        Some(true), // With prefiltering
    )
    .unwrap();

    let test_sequence = create_sequence_of_length(127);

    c.bench_function("single_sequence", |b| {
        b.iter(|| annotator.number_sequence(black_box(&test_sequence)))
    });
}

/// Benchmark algorithmic complexity by sequence length
fn benchmark_sequence_length_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("sequence_length_scaling");

    // Single chain for clean O(n²) analysis
    let annotator = Annotator::new(
        Scheme::IMGT,
        vec![Chain::IGH],
        None,
        Some(false), // No prefiltering for algorithm analysis
    )
    .unwrap();

    // Test genuinely different sequence lengths
    for length in [50, 100, 200, 300, 500, 750, 1000].iter() {
        let sequence = create_sequence_of_length(*length);

        group.throughput(Throughput::Elements(*length as u64));
        group.bench_with_input(
            BenchmarkId::new("needleman_wunsch", length),
            length,
            |b, _| b.iter(|| annotator.number_sequence(black_box(&sequence))),
        );
    }
    group.finish();
}

/// Benchmark batch processing with different sizes
fn benchmark_batch_processing(c: &mut Criterion) {
    let mut group = c.benchmark_group("batch_processing");
    group.measurement_time(Duration::from_secs(10)); // Longer measurement for accuracy

    let annotator = Annotator::new(
        Scheme::IMGT,
        vec![Chain::IGH, Chain::IGK, Chain::IGL],
        None,
        Some(true),
    )
    .unwrap();

    for batch_size in [10, 50, 100, 500, 1000].iter() {
        let sequences = create_realistic_test_sequences(*batch_size);

        group.throughput(Throughput::Elements(*batch_size as u64));

        // Sequential processing
        group.bench_with_input(
            BenchmarkId::new("sequential", batch_size),
            &sequences,
            |b, sequences| b.iter(|| annotator.number_sequences(black_box(sequences), false)),
        );

        // Parallel processing
        group.bench_with_input(
            BenchmarkId::new("parallel", batch_size),
            &sequences,
            |b, sequences| b.iter(|| annotator.number_sequences(black_box(sequences), true)),
        );
    }
    group.finish();
}

/// Benchmark prefiltering impact
fn benchmark_prefiltering_impact(c: &mut Criterion) {
    let mut group = c.benchmark_group("prefiltering_impact");

    let sequences = create_realistic_test_sequences(100);

    // Test different chain configurations
    let chain_configs = [
        ("single_chain", vec![Chain::IGH]),
        ("two_chains", vec![Chain::IGH, Chain::IGK]),
        ("three_chains", vec![Chain::IGH, Chain::IGK, Chain::IGL]),
    ];

    for (config_name, chains) in chain_configs.iter() {
        // Without prefiltering
        let annotator_no_filter =
            Annotator::new(Scheme::IMGT, chains.clone(), None, Some(false)).unwrap();

        group.bench_with_input(
            BenchmarkId::new("no_prefilter", config_name),
            &sequences,
            |b, sequences| {
                b.iter(|| annotator_no_filter.number_sequences(black_box(sequences), true))
            },
        );

        // With prefiltering
        let annotator_with_filter =
            Annotator::new(Scheme::IMGT, chains.clone(), None, Some(true)).unwrap();

        group.bench_with_input(
            BenchmarkId::new("with_prefilter", config_name),
            &sequences,
            |b, sequences| {
                b.iter(|| annotator_with_filter.number_sequences(black_box(sequences), true))
            },
        );
    }
    group.finish();
}

/// Benchmark threading efficiency at different work sizes
fn benchmark_threading_efficiency(c: &mut Criterion) {
    let mut group = c.benchmark_group("threading_efficiency");

    let annotator = Annotator::new(
        Scheme::IMGT,
        vec![Chain::IGH, Chain::IGK, Chain::IGL],
        None,
        Some(true),
    )
    .unwrap();

    // Test different amounts of work to find parallel break-even point
    for work_size in [1, 2, 5, 10, 20, 50, 100].iter() {
        let sequences = create_realistic_test_sequences(*work_size);

        group.throughput(Throughput::Elements(*work_size as u64));

        group.bench_with_input(
            BenchmarkId::new("sequential", work_size),
            &sequences,
            |b, sequences| b.iter(|| annotator.number_sequences(black_box(sequences), false)),
        );

        group.bench_with_input(
            BenchmarkId::new("parallel", work_size),
            &sequences,
            |b, sequences| b.iter(|| annotator.number_sequences(black_box(sequences), true)),
        );
    }
    group.finish();
}

/// Benchmark core Needleman-Wunsch algorithm directly
fn benchmark_needleman_wunsch_core(c: &mut Criterion) {
    let mut group = c.benchmark_group("needleman_wunsch_core");

    // Create a single scheme for direct algorithm testing
    let annotator = Annotator::new(Scheme::IMGT, vec![Chain::IGH], None, Some(false)).unwrap();

    // Test with different sequence complexities
    let test_cases = [
        ("short", create_sequence_of_length(80)),
        ("medium", create_sequence_of_length(150)),
        ("long", create_sequence_of_length(300)),
        ("very_long", create_sequence_of_length(600)),
    ];

    for (name, sequence) in test_cases.iter() {
        group.throughput(Throughput::Elements(sequence.len() as u64));
        group.bench_with_input(
            BenchmarkId::new("core_algorithm", name),
            sequence,
            |b, seq| b.iter(|| annotator.number_sequence(black_box(seq))),
        );
    }
    group.finish();
}

/// Benchmark different numbering schemes (IMGT vs KABAT)
fn benchmark_scheme_comparison(c: &mut Criterion) {
    let mut group = c.benchmark_group("scheme_comparison");
    group.measurement_time(Duration::from_secs(10));

    let sequences = create_realistic_test_sequences(50);

    // IMGT annotator
    let imgt_annotator = Annotator::new(
        Scheme::IMGT,
        vec![Chain::IGH, Chain::IGK, Chain::IGL],
        None,
        Some(true),
    )
    .unwrap();

    // KABAT annotator
    let kabat_annotator = Annotator::new(
        Scheme::KABAT,
        vec![Chain::IGH, Chain::IGK, Chain::IGL],
        None,
        Some(true),
    )
    .unwrap();

    group.throughput(Throughput::Elements(sequences.len() as u64));

    group.bench_function("imgt_sequential", |b| {
        b.iter(|| imgt_annotator.number_sequences(black_box(&sequences), false))
    });

    group.bench_function("imgt_parallel", |b| {
        b.iter(|| imgt_annotator.number_sequences(black_box(&sequences), true))
    });

    group.bench_function("kabat_sequential", |b| {
        b.iter(|| kabat_annotator.number_sequences(black_box(&sequences), false))
    });

    group.bench_function("kabat_parallel", |b| {
        b.iter(|| kabat_annotator.number_sequences(black_box(&sequences), true))
    });

    group.finish();
}

/// Benchmark different chain configurations
fn benchmark_chain_configurations(c: &mut Criterion) {
    let mut group = c.benchmark_group("chain_configurations");

    let test_sequence = create_sequence_of_length(200);

    let configs = [
        ("heavy_only", vec![Chain::IGH]),
        ("heavy_kappa", vec![Chain::IGH, Chain::IGK]),
        ("heavy_lambda", vec![Chain::IGH, Chain::IGL]),
        ("all_chains", vec![Chain::IGH, Chain::IGK, Chain::IGL]),
    ];

    for (name, chains) in configs.iter() {
        let annotator = Annotator::new(
            Scheme::IMGT,
            chains.clone(),
            None,
            Some(false), // No prefiltering for clean comparison
        )
        .unwrap();

        group.bench_function(*name, |b| {
            b.iter(|| annotator.number_sequence(black_box(&test_sequence)))
        });
    }

    group.finish();
}

/// Benchmark memory scaling with sequence length
fn benchmark_memory_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("memory_scaling");

    // Test with sequences of different lengths to see memory scaling
    let lengths = [100, 200, 400, 800, 1200];

    for &length in &lengths {
        let sequence = create_sequence_of_length(length);

        let annotator = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGH],
            None,
            Some(false), // No prefiltering for clean memory analysis
        )
        .unwrap();

        group.throughput(Throughput::Elements(length as u64));

        group.bench_with_input(
            BenchmarkId::new("sequence_length", length),
            &sequence,
            |b, seq| b.iter(|| annotator.number_sequence(black_box(seq))),
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    benchmark_single_sequence,
    benchmark_sequence_length_scaling,
    benchmark_needleman_wunsch_core,
    benchmark_prefiltering_impact,
    benchmark_threading_efficiency,
    benchmark_batch_processing,
    benchmark_scheme_comparison,
    benchmark_chain_configurations,
    benchmark_memory_scaling,
);
criterion_main!(benches);
