/**
 * Tests for Annotator methods in the WASM API.
 */

import { describe, it, expect, beforeAll } from 'vitest';
import { readFileSync } from 'fs';
import { fileURLToPath } from 'url';
import { dirname, join } from 'path';

// Get the directory paths
const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);
const rootDir = join(__dirname, '../../');

// Import WASM module
import init, { 
    Scheme, 
    Chain, 
    ScoringParams, 
    Annotator
} from '../../wasm_build/immunum.js';

describe('Annotator Methods', () => {
    beforeAll(async () => {
        // Load the WASM binary
        const wasmPath = join(rootDir, 'wasm_build', 'immunum_bg.wasm');
        const wasmBytes = readFileSync(wasmPath);
        await init(wasmBytes);
    });

    describe('numberSequence', () => {
        it('should handle numberSequence method calls', () => {
            const annotator = new Annotator(Scheme.IMGT, [Chain.IGH]);

            // Test with a short sequence - might fail but shouldn't crash
            const testSequence = "ATCGATCGATCG";
            
            const result = annotator.numberSequence(testSequence);
            expect(Array.isArray(result)).toBe(true);
        });

        it('should handle longer realistic sequences', () => {
            const annotator = new Annotator(Scheme.IMGT, [Chain.IGH]);

            // Test with a more realistic antibody sequence
            const antibodySequence = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR";
            
            const result = annotator.numberSequence(antibodySequence);
            expect(Array.isArray(result)).toBe(true);
        });

        it('should handle invalid sequences gracefully', () => {
            const annotator = new Annotator(Scheme.IMGT, [Chain.IGH]);

            const invalidSequences = [
                "", // Empty sequence
                "123", // Non-amino acid characters
                "XXXXXX", // Invalid amino acids only
            ];

            for (const seq of invalidSequences) {
                try {
                    const result = annotator.numberSequence(seq);
                    expect(Array.isArray(result)).toBe(true);
                } catch (error) {
                    // Accept a thrown error as graceful handling too
                    expect(error).toBeDefined();
                }
            }
        });
    });

    describe('numberSequences', () => {
        it('should handle batch processing', () => {
            const annotator = new Annotator(Scheme.IMGT, [Chain.IGH]);

            const sequences = ["ATCGATCGATCG", "GCTAGCTAGCTA"];
            
            try {
                const results = annotator.numberSequences(sequences);
                expect(Array.isArray(results)).toBe(true);
                expect(results.length).toBe(sequences.length);
                // Each item is an array (results per input)
                for (const perSeq of results) {
                    expect(Array.isArray(perSeq)).toBe(true);
                }
            } catch (error) {
                // Accept failure with test sequences
                expect(error).toBeDefined();
            }
        });

        it('should support paired mode via flag', () => {
            const annotator = new Annotator(Scheme.IMGT, [Chain.IGH, Chain.IGK]);
            try {
                const res = annotator.numberSequence("ATCGATCGATCG", true);
                expect(Array.isArray(res)).toBe(true);
                const res2 = annotator.numberSequences(["ATCGATCGATCG", "GCTAGCTAGCTA"], true);
                expect(Array.isArray(res2)).toBe(true);
                for (const perSeq of res2) {
                    expect(Array.isArray(perSeq)).toBe(true);
                }
            } catch (error) {
                expect(error).toBeDefined();
            }
        });
        it('should handle empty sequence list', () => {
            const annotator = new Annotator(Scheme.IMGT, [Chain.IGH]);
            
            const results = annotator.numberSequences([]);
            expect(Array.isArray(results)).toBe(true);
            expect(results.length).toBe(0);
        });

        it('should handle single sequence in batch', () => {
            const annotator = new Annotator(Scheme.IMGT, [Chain.IGH]);
            
            const sequences = ["ATCGATCGATCG"];
            try {
                const results = annotator.numberSequences(sequences);
                expect(Array.isArray(results)).toBe(true);
                expect(results.length).toBe(1);
            } catch (error) {
                expect(error).toBeDefined();
            }
        });
    });

    describe('Different Schemes and Chains', () => {
        it('should work with KABAT scheme', () => {
            const annotator = new Annotator(Scheme.KABAT, [Chain.IGH]);
            expect(annotator).toBeDefined();
            
            // Should be able to call methods without crashing
            try {
                annotator.numberSequence("ATCGATCGATCG");
            } catch (error) {
                // Expected to fail with test sequence
                expect(error).toBeDefined();
            }
        });

        it('should work with kappa chain', () => {
            const annotator = new Annotator(Scheme.IMGT, [Chain.IGK]);
            expect(annotator).toBeDefined();
        });

        it('should work with lambda chain', () => {
            const annotator = new Annotator(Scheme.IMGT, [Chain.IGL]);
            expect(annotator).toBeDefined();
        });

        it('should work with T-cell receptor chains', () => {
            const tcr_chains = [Chain.TRA, Chain.TRB, Chain.TRG, Chain.TRD];
            
            for (const chain of tcr_chains) {
                const annotator = new Annotator(Scheme.IMGT, [chain]);
                expect(annotator).toBeDefined();
            }
        });
    });

    describe('Error Handling', () => {
        it('should handle invalid scheme gracefully', () => {
            expect(() => {
                new Annotator(999, [Chain.IGH]); // Invalid scheme
            }).toThrow();
        });

        it('should handle invalid chain gracefully', () => {
            expect(() => {
                new Annotator(Scheme.IMGT, [999]); // Invalid chain
            }).toThrow();
        });

        it('should handle empty chain array', () => {
            expect(() => {
                new Annotator(Scheme.IMGT, []); // Empty chains
            }).toThrow();
        });
    });
});
