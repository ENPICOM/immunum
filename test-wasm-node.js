#!/usr/bin/env node

import { readFileSync } from 'fs';
import { fileURLToPath } from 'url';
import { dirname, join } from 'path';

// Get the directory of this script
const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);

// Import the WASM module
import init, { numberSequence } from './wasm_build/immunum.js';

async function testWasmInNode() {
    try {
        console.log('🧬 Testing Immunum WASM binary in Node.js...\n');
        
        // Load the WASM binary directly from file system
        const wasmPath = join(__dirname, 'wasm_build', 'immunum_bg.wasm');
        const wasmBytes = readFileSync(wasmPath);
        
        // Initialize the WASM module with the binary data
        await init(wasmBytes);
        
        console.log('✅ WASM module initialized successfully!\n');
        
        // Test cases
        const testCases = [
            {
                name: "Heavy chain (IGH) with IMGT scheme",
                sequence: "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR",
                scheme: "imgt",
                chains: ["igh"]
            },
            {
                name: "Kappa chain (IGK) with Kabat scheme", 
                sequence: "DIQMTQSPSSLSASVGDRVTITCRASQGIRNYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQRYNRAPYTFGQGTKVEIK",
                scheme: "kabat",
                chains: ["igk"]
            },
            {
                name: "Multiple chains test",
                sequence: "ACDEFGHIKLMNPQRSTVWY",
                scheme: "imgt", 
                chains: ["igh", "igk"]
            }
        ];
        
        // Run test cases
        for (let i = 0; i < testCases.length; i++) {
            const test = testCases[i];
            console.log(`🧪 Test ${i + 1}: ${test.name}`);
            console.log(`   Sequence: ${test.sequence.substring(0, 50)}${test.sequence.length > 50 ? '...' : ''}`);
            console.log(`   Scheme: ${test.scheme}`);
            console.log(`   Chains: ${test.chains.join(', ')}`);
            
            try {
                const result = numberSequence(test.sequence, test.scheme, test.chains);
                console.log(`   ✅ Result:\n${result}\n`);
            } catch (error) {
                console.log(`   ❌ Error: ${error.message}\n`);
            }
        }
        
        // Test error handling
        console.log('🧪 Testing error handling...');
        try {
            numberSequence("INVALID", "invalid_scheme", ["invalid_chain"]);
        } catch (error) {
            console.log(`   ✅ Error handling works: ${error.message}\n`);
        }
        
        console.log('🎉 All tests completed!');
        
    } catch (error) {
        console.error('❌ Failed to initialize WASM module:', error);
        process.exit(1);
    }
}

// Run the tests
testWasmInNode();
