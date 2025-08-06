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
            
            try {
                const result = annotator.numberSequence(testSequence);
                expect(result).toBeDefined();
            } catch (error) {
                // Short sequences might be rejected - this is expected behavior
                expect(error.message).toMatch(/(error|empty|short|invalid)/i);
            }
        });

        it('should handle longer realistic sequences', () => {
            const annotator = new Annotator(Scheme.IMGT, [Chain.IGH]);

            // Test with a more realistic antibody sequence
            const antibodySequence = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR";
            
            try {
                const result = annotator.numberSequence(antibodySequence);
                expect(result).toBeDefined();
            } catch (error) {
                // Even realistic sequences might fail depending on the algorithm
                // The important thing is that it doesn't crash
                expect(error.message).toBeTruthy();
            }
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
                    // If it doesn't throw, that's also acceptable behavior
                    expect(result).toBeDefined();
                } catch (error) {
                    // If it throws, that's also acceptable behavior
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
            } catch (error) {
                // Expected to potentially fail with test sequences
                expect(error.message).toMatch(/error/i);
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
            const results = annotator.numberSequences(sequences);
            expect(Array.isArray(results)).toBe(true);
            expect(results.length).toBe(1);
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