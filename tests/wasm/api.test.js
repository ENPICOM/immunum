/**
 * Basic tests for the immunum WASM API.
 * 
 * This module contains tests for the WASM bindings functionality,
 * mirroring the Python test structure.
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
    createAnnotator,
    Scheme,
    Chain,
    Annotator
} from '../../wasm_build/immunum.js';

describe('Immunum WASM API', () => {
    beforeAll(async () => {
        // Load the WASM binary
        const wasmPath = join(rootDir, 'wasm_build', 'immunum_bg.wasm');
        const wasmBytes = readFileSync(wasmPath);
        await init(wasmBytes);
    });

    describe('Enums', () => {
        it('should have accessible scheme enums', () => {
            expect(Scheme.IMGT).toBeDefined();
            expect(Scheme.KABAT).toBeDefined();
            expect(typeof Scheme.IMGT).toBe('number');
            expect(typeof Scheme.KABAT).toBe('number');
        });

        it('should have accessible chain enums', () => {
            expect(Chain.IGH).toBeDefined();
            expect(Chain.IGK).toBeDefined();
            expect(Chain.IGL).toBeDefined();
            expect(Chain.TRA).toBeDefined();
            expect(Chain.TRB).toBeDefined();
            expect(Chain.TRG).toBeDefined();
            expect(Chain.TRD).toBeDefined();
            
            // Verify they are numbers (WASM enums are numeric)
            expect(typeof Chain.IGH).toBe('number');
            expect(typeof Chain.IGK).toBe('number');
        });
    });

    describe('Sequence Numbering', () => {
        let annotator;

        beforeAll(() => {
            annotator = new Annotator(Scheme.IMGT, [Chain.IGH]);
        });

        it('should create numbering results from sequence', async () => {
            const heavyChainSequence = 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCAR';

            const result = annotator.numberSequence(heavyChainSequence, 1);
            expect(result).toBeDefined();
            expect(typeof result).toBe('string');

            // Should be valid JSON
            try {
                const parsed = JSON.parse(result);
                expect(Array.isArray(parsed)).toBe(true);

                if (parsed.length > 0) {
                    const chainResult = parsed[0];
                    expect(chainResult.numbers).toBeDefined();
                    expect(chainResult.confidence).toBeDefined();
                    expect(chainResult.chain).toBeDefined();
                    expect(chainResult.scheme).toBeDefined();
                    expect(chainResult.start).toBeDefined();
                    expect(chainResult.end).toBeDefined();

                    expect(Array.isArray(chainResult.numbers)).toBe(true);
                    expect(typeof chainResult.confidence).toBe('number');
                    expect(chainResult.confidence).toBeGreaterThan(0);
                }
            } catch (e) {
                // Check if it's an error result
                expect(result).toMatch(/^{"error":/);
            }
        });
    });

    describe('Annotator', () => {
        it('should create Annotator successfully', () => {
            const annotator = new Annotator(Scheme.IMGT, [Chain.IGH]);
            expect(annotator).toBeDefined();
        });

        it('should create Annotator with prefiltering enabled by default', () => {
            const annotator = new Annotator(
                Scheme.IMGT,
                [Chain.IGH, Chain.IGK, Chain.IGL]
            );
            expect(annotator).toBeDefined();
        });

        it('should create Annotator with prefiltering explicitly disabled', () => {
            const annotator = new Annotator(
                Scheme.IMGT,
                [Chain.IGH, Chain.IGK, Chain.IGL],
                true, // disable_prefiltering = true
                0.7   // min_confidence
            );
            expect(annotator).toBeDefined();
        });

        it('should create Annotator with custom confidence threshold', () => {
            const annotator = new Annotator(
                Scheme.IMGT,
                [Chain.IGH],
                false, // disable_prefiltering = false (prefiltering enabled)
                0.8    // min_confidence
            );
            expect(annotator).toBeDefined();
        });

        it('should create Annotator with multiple chains', () => {
            const annotator = new Annotator(
                Scheme.KABAT,
                [Chain.IGH, Chain.IGK]
            );
            expect(annotator).toBeDefined();
        });

        it('should number a single sequence', () => {
            const annotator = new Annotator(Scheme.IMGT, [Chain.IGH]);
            const sequence = 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCAR';

            const result = annotator.numberSequence(sequence, 1);
            expect(result).toBeDefined();
            expect(typeof result).toBe('string');

            // Should be valid JSON
            try {
                const parsed = JSON.parse(result);
                expect(Array.isArray(parsed)).toBe(true);
            } catch (e) {
                // Check if it's an error result
                expect(result).toMatch(/^{"error":/);
            }
        });

        it('should number multiple sequences', () => {
            const annotator = new Annotator(Scheme.IMGT, [Chain.IGH, Chain.IGK]);
            const sequences = [
                'QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCAR',
                'DIVMTQSPDSLAVSLGERATINCKASQSVTNDVAWYQQKPGQPPKLLIYYASNRYTGVPDRFSGSGSGTDFTLTISSLQAEDVAVYYCQQDYSSPYT'
            ];

            // New API takes JSON string input
            const sequencesJson = JSON.stringify(sequences);
            const result = annotator.numberSequences(sequencesJson, 2);
            expect(result).toBeDefined();
            expect(typeof result).toBe('string');

            // Should return a JSON array of results
            try {
                const parsed = JSON.parse(result);
                expect(Array.isArray(parsed)).toBe(true);
                expect(parsed.length).toBe(2);
            } catch (e) {
                // Check if it's an error result
                expect(result).toMatch(/^{"error":/);
            }
        });
    });

    describe('API Structure', () => {
        it('should have all expected exports available', () => {
            // Test that all main API components are available
            expect(createAnnotator).toBeDefined();
            expect(Scheme).toBeDefined();
            expect(Chain).toBeDefined();
            expect(Annotator).toBeDefined();

            // Test that they are the right types
            expect(typeof createAnnotator).toBe('function');
            expect(typeof Annotator).toBe('function'); // Constructor
        });

        it('should create annotator using createAnnotator helper', () => {
            const annotator = createAnnotator();
            expect(annotator).toBeDefined();
        });

        it('should create annotator with custom parameters using createAnnotator', () => {
            const annotator = createAnnotator(
                Scheme.KABAT,
                [Chain.IGH, Chain.IGK],
                true, // disable prefiltering
                0.8   // min confidence
            );
            expect(annotator).toBeDefined();
        });
    });
});