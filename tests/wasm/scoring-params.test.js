/**
 * Comprehensive tests for ScoringParams in the WASM API.
 * 
 * These tests verify that the Rust ScoringParams struct works correctly
 * through WASM bindings with automatic getter/setter generation.
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
    ScoringParams 
} from '../../wasm_build/immunum.js';

describe('ScoringParams WASM Integration', () => {
    beforeAll(async () => {
        // Load the WASM binary
        const wasmPath = join(rootDir, 'wasm_build', 'immunum_bg.wasm');
        const wasmBytes = readFileSync(wasmPath);
        await init(wasmBytes);
    });

    describe('Constructor and Default Values', () => {
        it('should create with default values', () => {
            const params = new ScoringParams();
            
            // Test that all fields exist and have reasonable default values
            expect(params.gap_pen_cp).toBe(55.0);
            expect(params.gap_pen_fr).toBe(26.0);
            expect(params.gap_pen_ip).toBe(1.5);
            expect(params.gap_pen_op).toBe(1.0);
            expect(params.gap_pen_cdr).toBe(2.5);
            expect(params.gap_pen_other).toBe(11.0);
            expect(params.cdr_increase).toBe(0.5);
            expect(params.pen_leap_insertion_point_imgt).toBe(1.0);
            expect(params.pen_leap_insertion_point_kabat).toBe(10.0);
        });

        it('should create with partial custom values', () => {
            const params = new ScoringParams(100.0, 50.0);
            
            expect(params.gap_pen_cp).toBe(100.0);
            expect(params.gap_pen_fr).toBe(50.0);
            // Other values should be defaults
            expect(params.gap_pen_ip).toBe(1.5);
            expect(params.gap_pen_op).toBe(1.0);
        });

        it('should create with all custom values', () => {
            const params = new ScoringParams(
                100.0, // gap_pen_cp
                200.0, // gap_pen_fr  
                3.0,   // gap_pen_ip
                4.0,   // gap_pen_op
                5.0,   // gap_pen_cdr
                6.0,   // gap_pen_other
                7.0,   // cdr_increase
                8.0,   // pen_leap_insertion_point_imgt
                9.0    // pen_leap_insertion_point_kabat
            );
            
            expect(params.gap_pen_cp).toBe(100.0);
            expect(params.gap_pen_fr).toBe(200.0);
            expect(params.gap_pen_ip).toBe(3.0);
            expect(params.gap_pen_op).toBe(4.0);
            expect(params.gap_pen_cdr).toBe(5.0);
            expect(params.gap_pen_other).toBe(6.0);
            expect(params.cdr_increase).toBe(7.0);
            expect(params.pen_leap_insertion_point_imgt).toBe(8.0);
            expect(params.pen_leap_insertion_point_kabat).toBe(9.0);
        });
    });

    describe('Default Scoring Parameters Function', () => {
        it('should return ScoringParams with default values', () => {
            const params = defaultScoringParams();
            
            expect(params).toBeInstanceOf(ScoringParams);
            expect(params.gap_pen_cp).toBe(55.0);
            expect(params.gap_pen_fr).toBe(26.0);
        });

        it('should return equivalent values to constructor defaults', () => {
            const defaultParams = defaultScoringParams();
            const constructorParams = new ScoringParams();
            
            expect(defaultParams.gap_pen_cp).toBe(constructorParams.gap_pen_cp);
            expect(defaultParams.gap_pen_fr).toBe(constructorParams.gap_pen_fr);
            expect(defaultParams.gap_pen_ip).toBe(constructorParams.gap_pen_ip);
            expect(defaultParams.gap_pen_op).toBe(constructorParams.gap_pen_op);
            expect(defaultParams.gap_pen_cdr).toBe(constructorParams.gap_pen_cdr);
            expect(defaultParams.gap_pen_other).toBe(constructorParams.gap_pen_other);
            expect(defaultParams.cdr_increase).toBe(constructorParams.cdr_increase);
            expect(defaultParams.pen_leap_insertion_point_imgt).toBe(constructorParams.pen_leap_insertion_point_imgt);
            expect(defaultParams.pen_leap_insertion_point_kabat).toBe(constructorParams.pen_leap_insertion_point_kabat);
        });
    });

    describe('Getters and Setters', () => {
        it('should have working getters for all fields', () => {
            const params = new ScoringParams();
            
            // All getters should return numbers
            expect(typeof params.gap_pen_cp).toBe('number');
            expect(typeof params.gap_pen_fr).toBe('number');
            expect(typeof params.gap_pen_ip).toBe('number');
            expect(typeof params.gap_pen_op).toBe('number');
            expect(typeof params.gap_pen_cdr).toBe('number');
            expect(typeof params.gap_pen_other).toBe('number');
            expect(typeof params.cdr_increase).toBe('number');
            expect(typeof params.pen_leap_insertion_point_imgt).toBe('number');
            expect(typeof params.pen_leap_insertion_point_kabat).toBe('number');
        });

        it('should have working setters for all fields', () => {
            const params = new ScoringParams();
            
            // Test each setter
            params.gap_pen_cp = 999.1;
            expect(params.gap_pen_cp).toBe(999.1);
            
            params.gap_pen_fr = 999.2;
            expect(params.gap_pen_fr).toBe(999.2);
            
            params.gap_pen_ip = 999.3;
            expect(params.gap_pen_ip).toBe(999.3);
            
            params.gap_pen_op = 999.4;
            expect(params.gap_pen_op).toBe(999.4);
            
            params.gap_pen_cdr = 999.5;
            expect(params.gap_pen_cdr).toBe(999.5);
            
            params.gap_pen_other = 999.6;
            expect(params.gap_pen_other).toBe(999.6);
            
            params.cdr_increase = 999.7;
            expect(params.cdr_increase).toBe(999.7);
            
            params.pen_leap_insertion_point_imgt = 999.8;
            expect(params.pen_leap_insertion_point_imgt).toBe(999.8);
            
            params.pen_leap_insertion_point_kabat = 999.9;
            expect(params.pen_leap_insertion_point_kabat).toBe(999.9);
        });

        it('should handle negative values', () => {
            const params = new ScoringParams();
            
            params.gap_pen_cp = -10.5;
            expect(params.gap_pen_cp).toBe(-10.5);
        });

        it('should handle zero values', () => {
            const params = new ScoringParams();
            
            params.gap_pen_cp = 0.0;
            expect(params.gap_pen_cp).toBe(0.0);
        });

        it('should handle very large values', () => {
            const params = new ScoringParams();
            
            params.gap_pen_cp = 1e10;
            expect(params.gap_pen_cp).toBe(1e10);
        });

        it('should handle very small values', () => {
            const params = new ScoringParams();
            
            params.gap_pen_cp = 1e-10;
            expect(params.gap_pen_cp).toBe(1e-10);
        });
    });

    describe('Memory Management', () => {
        it('should be able to create many instances', () => {
            const instances = [];
            
            // Create many instances to test memory handling
            for (let i = 0; i < 100; i++) {
                instances.push(new ScoringParams(i * 1.0));
            }
            
            // Verify they all work
            for (let i = 0; i < 100; i++) {
                expect(instances[i].gap_pen_cp).toBe(i * 1.0);
            }
            
            // Clean up (call free() if available)
            instances.forEach(instance => {
                if (typeof instance.free === 'function') {
                    instance.free();
                }
            });
        });

        it('should maintain independent state across instances', () => {
            const params1 = new ScoringParams();
            const params2 = new ScoringParams();
            
            // Modify one instance
            params1.gap_pen_cp = 123.45;
            
            // Other instance should be unchanged
            expect(params2.gap_pen_cp).toBe(55.0); // Default value
            expect(params1.gap_pen_cp).toBe(123.45);
        });
    });

    describe('Integration with Other Components', () => {
        it('should work as parameter to Annotator', async () => {
            // Import additional dependencies
            const { Annotator, Scheme, Chain } = await import('../../wasm_build/immunum.js');
            
            const customParams = new ScoringParams(99.0, 88.0);
            
            // Should not throw when passed to Annotator
            expect(() => {
                const annotator = new Annotator(Scheme.IMGT, [Chain.IGH], customParams);
                expect(annotator).toBeDefined();
            }).not.toThrow();
        });
    });
});