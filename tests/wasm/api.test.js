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
    defaultScoringParams, 
    Scheme, 
    Chain, 
    ScoringParams, 
    Annotator,
    AnnotationResult 
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

    describe('ScoringParams', () => {
        it('should create default ScoringParams', () => {
            const params = new ScoringParams();
            expect(params).toBeDefined();
            expect(params.gap_pen_cp).toBeGreaterThan(0);
            expect(typeof params.gap_pen_cp).toBe('number');
        });

        it('should create ScoringParams with custom values', () => {
            const customParams = new ScoringParams(60.0);
            expect(customParams).toBeDefined();
            expect(customParams.gap_pen_cp).toBe(60.0);
        });

        it('should have working getters and setters', () => {
            const params = new ScoringParams();
            const originalValue = params.gap_pen_cp;
            
            // Test setter
            params.gap_pen_cp = 123.45;
            expect(params.gap_pen_cp).toBe(123.45);
            
            // Test that it actually changed
            expect(params.gap_pen_cp).not.toBe(originalValue);
        });

        it('should have all expected properties', () => {
            const params = new ScoringParams();
            
            expect(params.gap_pen_cp).toBeDefined();
            expect(params.gap_pen_fr).toBeDefined();
            expect(params.gap_pen_ip).toBeDefined();
            expect(params.gap_pen_op).toBeDefined();
            expect(params.gap_pen_cdr).toBeDefined();
            expect(params.gap_pen_other).toBeDefined();
            expect(params.cdr_increase).toBeDefined();
            expect(params.pen_leap_insertion_point_imgt).toBeDefined();
            expect(params.pen_leap_insertion_point_kabat).toBeDefined();
        });

        it('should create default scoring parameters via function', () => {
            const params = defaultScoringParams();
            expect(params).toBeDefined();
            expect(params.gap_pen_cp).toBeGreaterThan(0);
            expect(typeof params.gap_pen_cp).toBe('number');
        });
    });

    describe('Annotator', () => {
        it('should create Annotator successfully', () => {
            const annotator = new Annotator(Scheme.IMGT, [Chain.IGH]);
            expect(annotator).toBeDefined();
        });

        it('should create Annotator with custom scoring parameters', () => {
            const customParams = new ScoringParams(60.0);
            const annotator = new Annotator(
                Scheme.IMGT, 
                [Chain.IGH], 
                customParams
            );
            expect(annotator).toBeDefined();
        });

        it('should create Annotator with prefiltering enabled', () => {
            const annotator = new Annotator(
                Scheme.IMGT, 
                [Chain.IGH, Chain.IGK, Chain.IGL], 
                null, 
                true
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
    });

    describe('API Structure', () => {
        it('should have all expected exports available', () => {
            // Test that all main API components are available
            expect(defaultScoringParams).toBeDefined();
            expect(Scheme).toBeDefined();
            expect(Chain).toBeDefined();
            expect(ScoringParams).toBeDefined();
            expect(Annotator).toBeDefined();
            expect(AnnotationResult).toBeDefined();
            
            // Test that they are the right types
            expect(typeof defaultScoringParams).toBe('function');
            expect(typeof ScoringParams).toBe('function'); // Constructor
            expect(typeof Annotator).toBe('function'); // Constructor
            expect(typeof AnnotationResult).toBe('function'); // Constructor
        });
    });
});