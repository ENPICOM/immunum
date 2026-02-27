use immunum::{Annotator, Chain};

const ALL_CHAINS: &[Chain] = &[
    Chain::IGH,
    Chain::IGK,
    Chain::IGL,
    Chain::TRA,
    Chain::TRB,
    Chain::TRG,
    Chain::TRD,
];

fn main() {
    println!("immunum - Antibody and TCR Numbering Tool\n");

    // Create annotator with all chain types for auto-detection
    let annotator = match Annotator::new(ALL_CHAINS, immunum::Scheme::IMGT) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("Error creating annotator: {}", e);
            return;
        }
    };

    println!("Loaded annotator with {} chain types\n", ALL_CHAINS.len());

    // Test sequences
    let test_sequences = vec![
        (
            "IGH",
            "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNAKN",
        ),
        (
            "IGK",
            "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLTFGGGTKVEIK",
        ),
    ];

    for (expected, sequence) in test_sequences {
        println!("Sequence: {}... (len={})", &sequence[..40], sequence.len());
        println!("Expected chain: {}\n", expected);

        match annotator.annotate(sequence) {
            Ok(result) => {
                println!("Detected chain: {}", result.chain);
                println!("Confidence: {:.2}", result.confidence);
                println!("Alignment score: {:.2}", result.alignment.score);

                let numbering = result.numbering(immunum::Scheme::IMGT);
                println!("\nNumbering (first 10 positions):");
                for (i, (aa, pos)) in sequence.chars().zip(numbering.iter()).take(10).enumerate() {
                    println!("  {}: {} -> {}", i + 1, aa, pos);
                }
                println!("\n{}\n", "=".repeat(60));
            }
            Err(e) => {
                eprintln!("Error annotating sequence: {}\n", e);
            }
        }
    }
}
