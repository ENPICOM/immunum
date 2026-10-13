---
title: Demo
hide:
  - navigation
  - toc
  # Previous/next would lead out of the tab, to the neighbouring tabs' pages
  - footer
---

<div class="immunum-hero" markdown>
  <img src="../assets/immunum_logo.svg" alt="immunum" class="immunum-hero-logo" />
  <div class="immunum-hero-text">
    <h1>immunum</h1>
    <p>Fast antibody and T-cell receptor numbering, in your browser.</p>
  </div>
</div>

## Try it in your browser

<div class="immunum-tool" id="immunum-tool">
  <form id="immunum-form" class="immunum-form">
    <label class="immunum-field immunum-field-seq">
      <span>Amino-acid sequence or FASTA</span>
      <textarea id="seq-input" rows="4" spellcheck="false" autocomplete="off" aria-describedby="seq-hint"
        placeholder="Paste a heavy, light or TCR chain, or an scFv with several domains..."></textarea>
      <small id="seq-hint" class="immunum-hint">Ctrl+Enter (⌘+Enter on a Mac) numbers the sequence. Of a FASTA paste, only the first record is numbered.</small>
    </label>

    <div class="immunum-examples">
      <span>Load example:</span>
      <button type="button" data-example="IGH">IGH</button>
      <button type="button" data-example="IGL">IGL</button>
      <button type="button" data-example="TRB">TRB</button>
      <button type="button" data-example="scFv">scFv (VH–linker–VL)</button>
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
          <label class="immunum-scheme-radio">
            <input type="radio" name="scheme" value="chothia">
            <span>Chothia</span>
          </label>
          <label class="immunum-scheme-radio">
            <input type="radio" name="scheme" value="martin">
            <span>Martin</span>
          </label>
          <label class="immunum-scheme-radio">
            <input type="radio" name="scheme" value="aho">
            <span>AHo</span>
          </label>
        </div>
      </fieldset>

      <fieldset class="immunum-field immunum-confidence-field">
        <legend>Min. confidence <span id="confidence-value">0.50</span></legend>
        <input type="range" id="min-confidence" min="0" max="1" step="0.01" value="0.5" autocomplete="off">
      </fieldset>

      <div class="immunum-run-wrap">
        <button type="submit" id="run-btn" class="immunum-run" disabled>
          Loading WASM…
        </button>
      </div>
    </div>
  </form>

  <div id="immunum-error" class="immunum-error" role="alert"></div>
  <p id="immunum-notice" class="immunum-notice" role="status"></p>

  <div id="immunum-result" class="immunum-result" hidden>
    <div class="immunum-result-header">
      <span class="immunum-badge" id="result-domain-count">—</span>
      <span class="immunum-badge" id="result-scheme">—</span>
    </div>

    <div class="immunum-actions">
      <button type="button" id="copy-tsv">Copy TSV</button>
      <button type="button" id="download-tsv">Download TSV</button>
      <button type="button" id="download-json">Download JSON</button>
      <button type="button" id="copy-link">Copy link</button>
      <span id="immunum-action-status" class="immunum-action-status" role="status"></span>
    </div>

    <div class="immunum-legend">
      <span class="region-prefix">Flank</span>
      <span class="region-fr1">FR1</span>
      <span class="region-cdr1">CDR1</span>
      <span class="region-fr2">FR2</span>
      <span class="region-cdr2">CDR2</span>
      <span class="region-fr3">FR3</span>
      <span class="region-cdr3">CDR3</span>
      <span class="region-fr4">FR4</span>
    </div>

    <div class="immunum-aligned">
      <span class="immunum-aligned-label">Sequence</span>
      <pre id="result-sequence"></pre>
    </div>

    <div id="result-domains"></div>
  </div>
  <div class="immunum-version" id="immunum-version" hidden></div>
</div>

## How to read the result

- **Sequence**: your whole input, coloured by the region each residue falls in: FR1 to FR4 and
  CDR1 to CDR3. Hover over a stretch to see its domain and region.
- **Flank**: grey residues outside every domain. A domain's prefix holds the residues before it,
  back to the previous domain or the start of the sequence (a leader, or the linker of an scFv);
  the last domain's postfix holds the residues after it. See
  [Prefix and postfix](concepts.md#prefix-and-postfix).
- **Query**: the residues of your input that the domain's alignment covers, counted from 1 and
  inclusive. The APIs report the same range 0-based, as `queryStart`/`queryEnd` (`query_start`/
  `query_end` in Python and the CLI).
- **Confidence**: how well the domain aligns to the chain it was assigned, from 0 to 1. A domain
  below the minimum confidence isn't numbered. See [Confidence](concepts.md#confidence).
- **Numbering**: one box per numbered residue, with its position in the chosen scheme above it.
  A position with a letter, such as `111A`, is an insertion: a residue the scheme has no plain
  number for. In IMGT CDR3 insertions grow from the middle, so `111A` follows `111` and `112A`
  comes before `112`; Kabat, Chothia and Martin append them after the number (`100`, `100A`,
  `100B`). Positions missing from the run, such as IMGT `10` in the IGH example, are gaps: the
  sequence has no residue there.
- **Compare schemes**: per domain, where each scheme puts the FR/CDR boundaries for that chain and
  which residues each scheme calls CDRs. Only IMGT numbers TCR chains.
- **Copy link**: a URL holding the sequence and settings; opening it numbers the sequence again.
  Downloads use the columns and keys of the CLI's `immunum number --all-domains` output, and the
  JSON adds each domain's segments.

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

<script type="module" src="../assets/immunum-webtool.js"></script>
