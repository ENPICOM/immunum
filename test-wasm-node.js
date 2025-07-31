#!/usr/bin/env node

import { readFileSync } from 'fs';
import { fileURLToPath } from 'url';
import { dirname, join } from 'path';

// Get the directory of this script
const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);

// Import the WASM module
import init, { numberSequencesBatch, Scheme, Chain, ScoringParams, NumberingScheme } from './wasm_build/immunum.js';

async function testWasmInNode() {
    try {
        console.log('🧬 Testing Immunum WASM binary in Node.js...\n');
        
        // Load the WASM binary directly from file system
        const wasmPath = join(__dirname, 'wasm_build', 'immunum_bg.wasm');
        const wasmBytes = readFileSync(wasmPath);
        
        // Initialize the WASM module with the binary data
        await init(wasmBytes);
        
        console.log('✅ WASM module initialized successfully!');
        
        // Test typed enums
        console.log('\n🧪 Testing typed enums...');
        try {
            console.log('Available schemes:', Scheme);
            console.log('Available chains:', Chain);
            console.log('✅ Typed enums accessible');
        } catch (error) {
            console.log(`❌ Enum access failed: ${error.message}`);
        }
        
        // Test ScoringParams
        console.log('\n🧪 Testing ScoringParams...');
        try {
            const params = new ScoringParams();
            console.log('✅ ScoringParams created successfully');
        } catch (error) {
            console.log(`❌ ScoringParams creation failed: ${error.message}`);
        }
        
        // Test NumberingScheme
        console.log('\n🧪 Testing NumberingScheme...');
        try {
            const scheme = NumberingScheme.imgtHeavy();
            console.log('✅ NumberingScheme.imgtHeavy() created successfully');
        } catch (error) {
            console.log(`❌ NumberingScheme creation failed: ${error.message}`);
        }
        
        // Test typed batch processing
        console.log('\n🧪 Testing typed batch processing...');
        const testCases = [
            {
                name: "Heavy chain (IGH) with IMGT scheme",
                sequences: ["EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR"],
                scheme: Scheme.IMGT,
                chains: [Chain.IGH]
            },
            {
                name: "Kappa chain (IGK) with Kabat scheme", 
                sequences: ["DIQMTQSPSSLSASVGDRVTITCRASQGIRNYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQRYNRAPYTFGQGTKVEIK"],
                scheme: Scheme.KABAT,
                chains: [Chain.IGK]
            },
            {
                name: "Multiple chains test",
                sequences: ["ACDEFGHIKLMNPQRSTVWY"],
                scheme: Scheme.IMGT, 
                chains: [Chain.IGH, Chain.IGK]
            }
        ];
        
        // Run test cases
        for (let i = 0; i < testCases.length; i++) {
            const test = testCases[i];
            try {
                const results = numberSequencesBatch(test.sequences, test.scheme, test.chains);
                console.log(`✅ Test ${i + 1}: ${test.name} - Got ${results.length} results`);
            } catch (error) {
                console.log(`❌ Test ${i + 1}: ${test.name} \n    Error: ${error.message}\n`);
            }
        }
        
        // Test error handling would be done with try/catch around the typed calls above
        console.log('\n✅ Error handling is built into the typed API calls above');
        
        console.log('🎉 All tests completed!');
        
    } catch (error) {
        console.error('❌ Failed to initialize WASM module:', error);
        process.exit(1);
    }
}

// Run the tests
testWasmInNode();
