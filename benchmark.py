#!/usr/bin/env python3
"""
Immunum Python API Performance Benchmark

This script benchmarks parallel processing and prefiltering performance
in the Immunum Python API.
"""

import time
import immunum
from typing import List

# Sample sequences for testing (mix of heavy and light chains)
SAMPLE_SEQUENCES = [
    # Heavy chain sequences
    "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARWGGRGSYAMDYWGQGTLVTVSS",
    "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDGGYCSGGSCYFDYWGQGTLVTVSS",
    "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGGIIPIFGTANYAQKFQGRVTITADESTSTAYMELSSLRSEDTAVYYCARSHYGMDVWGQGTTVTVSS",
    # Light chain sequences (kappa)
    "DIQMTQSPSSLSASVGDRVTITCRASQGIRNDLGWYQQKPGKAPKRLIYDASSLESGVPSRFSGSGSGTDFTFTISSLQPEDIATYYCQQSYSTPWTFGQGTKVEIK",
    "DIVMTQSHKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKALIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELK",
    "DIQMTQSPSSLSASVGDRVTITCRASQSISSWLAWYQQKPGKAPKLLIYDASSLESGVPSRFSGSGSGTDFTFTISSLQPEDIATYYCQQYYSTPLTFGQGTKVEIK",
    # Light chain sequences (lambda)
    "QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSKRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCSSYAGSNNLVFGAGTKLTVL",
    "SYELTQPPSVSVSPGQTARITCSGDALPKKYAYWYQQKPGQSPVLVIYDSKRPSGIPERFSGSNSGNTATLISRVEAGDEADYYCQSADSSGTYVFGTGTKVTVL",
]


def create_test_sequences(count: int) -> List[str]:
    """Generate a list of test sequences by cycling through samples."""
    sequences = []
    for i in range(count):
        sequences.append(SAMPLE_SEQUENCES[i % len(SAMPLE_SEQUENCES)])
    return sequences


def quick_test():
    """Quick test with 2000 sequences comparing prefiltering and parallel processing."""
    print("🚀 Immunum Performance Test - 200 Sequences")
    print("=" * 50)

    # Test with 200 sequences
    sequences = create_test_sequences(2000)
    print(f"Testing with {len(sequences)} sequences...")
    print()

    # Test 1: Without prefiltering
    print("📊 WITHOUT Prefiltering:")
    print("-" * 30)
    annotator_no_filter = immunum.Annotator(
        scheme=immunum.Scheme.IMGT,
        chains=[immunum.Chain.IGH, immunum.Chain.IGK, immunum.Chain.IGL],
        use_prefiltering=False,
    )

    start = time.time()
    sequential_results_no_filter = annotator_no_filter.number_sequences(sequences)
    seq_time_no_filter = time.time() - start

    start = time.time()
    parallel_results_no_filter = annotator_no_filter.number_sequences(
        sequences, parallel=True
    )
    par_time_no_filter = time.time() - start

    print(
        f"  Sequential: {seq_time_no_filter:.3f}s ({len(sequential_results_no_filter)} results)"
    )
    print(
        f"  Parallel:   {par_time_no_filter:.3f}s ({len(parallel_results_no_filter)} results)"
    )
    print(f"  Speedup:    {seq_time_no_filter / par_time_no_filter:.2f}x faster")
    print()

    # Test 2: With prefiltering
    print("📊 WITH Prefiltering:")
    print("-" * 30)
    annotator_with_filter = immunum.Annotator(
        scheme=immunum.Scheme.IMGT,
        chains=[immunum.Chain.IGH, immunum.Chain.IGK, immunum.Chain.IGL],
        use_prefiltering=True,
    )

    start = time.time()
    sequential_results_with_filter = annotator_with_filter.number_sequences(
        sequences, parallel=False
    )
    seq_time_with_filter = time.time() - start

    start = time.time()
    parallel_results_with_filter = annotator_with_filter.number_sequences(
        sequences, parallel=True
    )
    par_time_with_filter = time.time() - start

    print(
        f"  Sequential: {seq_time_with_filter:.3f}s ({len(sequential_results_with_filter)} results)"
    )
    print(
        f"  Parallel:   {par_time_with_filter:.3f}s ({len(parallel_results_with_filter)} results)"
    )
    print(f"  Speedup:    {seq_time_with_filter / par_time_with_filter:.2f}x faster")
    print()

    # Comparison
    print("🔍 Comparison Summary:")
    print("-" * 30)
    seq_prefilter_speedup = seq_time_no_filter / seq_time_with_filter
    par_prefilter_speedup = par_time_no_filter / par_time_with_filter

    print(f"Prefiltering benefit (Sequential): {seq_prefilter_speedup:.2f}x faster")
    print(f"Prefiltering benefit (Parallel):   {par_prefilter_speedup:.2f}x faster")
    print()
    print(f"Best performance: Parallel with prefiltering = {par_time_with_filter:.3f}s")
    print(
        f"Overall speedup vs sequential no prefiltering: {seq_time_no_filter / par_time_with_filter:.2f}x faster"
    )


if __name__ == "__main__":
    quick_test()
