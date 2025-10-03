import { defineConfig } from 'vitest/config';

export default defineConfig({
  test: {
    globals: true,
    environment: 'node',
    include: ['tests/**/*.{test,spec}.{js,mjs,ts}'],
    testTimeout: 10000, // 10 seconds for WASM loading
  },
});