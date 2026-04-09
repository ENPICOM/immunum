---
title: immunum
hide:
  - navigation
  - toc
---

<div class="immunum-hero" markdown>
  <img src="assets/immunum_logo.svg" alt="immunum" class="immunum-hero-logo" />
  <div class="immunum-hero-text">
    <h1>immunum</h1>
    <p>Fast antibody and T-cell receptor numbering, in your browser.</p>
  </div>
</div>

## Try it in your browser

<div class="immunum-tool" id="immunum-tool">
  <form id="immunum-form" class="immunum-form">
    <label class="immunum-field immunum-field-seq">
      <span>Amino-acid sequence</span>
      <textarea id="seq-input" rows="4" spellcheck="false" autocomplete="off"
        placeholder="Paste a heavy, light, or TCR chain sequence..."></textarea>
    </label>

    <div class="immunum-examples">
      <span>Load example:</span>
      <button type="button" data-example="IGH">IGH</button>
      <button type="button" data-example="IGL">IGL</button>
      <button type="button" data-example="TRB">TRB</button>
    </div>

    <div class="immunum-controls">
      <fieldset class="immunum-field immunum-chain-field">
        <legend>Chains</legend>
        <div id="chain" class="immunum-chain-grid"></div>
      </fieldset>

      <fieldset class="immunum-field immunum-scheme-field">
        <legend>Scheme</legend>
        <div class="immunum-scheme-grid" id="scheme">
          <label class="immunum-scheme-radio">
            <input type="radio" name="scheme" value="imgt" checked>
            <span>IMGT</span>
          </label>
          <label class="immunum-scheme-radio">
            <input type="radio" name="scheme" value="kabat">
            <span>Kabat</span>
          </label>
        </div>
      </fieldset>

      <div class="immunum-run-wrap">
        <button type="submit" id="run-btn" class="immunum-run" disabled>
          Loading WASM…
        </button>
      </div>
    </div>
  </form>

  <div id="immunum-error" class="immunum-error" hidden></div>

  <div id="immunum-result" class="immunum-result" hidden>
    <div class="immunum-result-header">
      <span class="immunum-badge" id="result-chain">—</span>
      <span class="immunum-badge" id="result-scheme">—</span>
      <span class="immunum-badge" id="result-range">—</span>
      <span class="immunum-badge" id="result-confidence-value">—</span>
    </div>

    <div class="immunum-legend">
      <span class="region-fr1">FR1</span>
      <span class="region-cdr1">CDR1</span>
      <span class="region-fr2">FR2</span>
      <span class="region-cdr2">CDR2</span>
      <span class="region-fr3">FR3</span>
      <span class="region-cdr3">CDR3</span>
      <span class="region-fr4">FR4</span>
    </div>

    <div id="result-grid" class="immunum-grid"></div>

    <div class="immunum-aligned">
      <span class="immunum-aligned-label">Aligned region</span>
      <pre id="result-aligned"></pre>
    </div>
  </div>
</div>

## Resources

<div class="immunum-cards">
  <a class="immunum-card" href="https://crates.io/crates/immunum" target="_blank" rel="noopener">
    <h3>crates.io</h3>
    <p>The Rust crate on crates.io.</p>
  </a>
  <a class="immunum-card" href="https://docs.rs/immunum/" target="_blank" rel="noopener">
    <h3>docs.rs</h3>
    <p>API documentation for the Rust crate.</p>
  </a>
  <a class="immunum-card" href="https://pypi.org/project/immunum/" target="_blank" rel="noopener">
    <h3>PyPI</h3>
    <p>The Python package.</p>
  </a>
  <a class="immunum-card" href="https://www.npmjs.com/package/immunum" target="_blank" rel="noopener">
    <h3>npm</h3>
    <p>The WASM package for JavaScript and TypeScript.</p>
  </a>
</div>

<script type="module" src="assets/immunum-webtool.js"></script>
