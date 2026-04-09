// Interactive WASM numbering tool for the immunum docs homepage.
import init, { Annotator } from "./wasm/immunum.js";

const EXAMPLES = {
  IGH: "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS",
  IGL: "SALTQPPAVSGTPGQRVTISCSGSDIGRRSVNWYQQFPGTAPKLLIYSNDQRPSVVPDRFSGSKSGTSASLAISGLQSEDEAEYYCAAWDDSLAVFGGGTQLTVGQPKA",
  TRB: "GVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASRPGLAGGRPEQYFGPGTRLTVTE",
};

const ALL_CHAINS = [
  { value: "H", label: "Heavy (IGH)" },
  { value: "K", label: "Kappa (IGK)" },
  { value: "L", label: "Lambda (IGL)" },
  { value: "A", label: "Alpha (TRA)" },
  { value: "B", label: "Beta (TRB)" },
  { value: "G", label: "Gamma (TRG)" },
  { value: "D", label: "Delta (TRD)" },
];

// Kabat is only defined for antibody chains.
const KABAT_CHAINS = new Set(["H", "K", "L"]);

function isChainAllowed(chain, scheme) {
  return scheme === "kabat" ? KABAT_CHAINS.has(chain) : true;
}

const REGIONS = ["fr1", "cdr1", "fr2", "cdr2", "fr3", "cdr3", "fr4"];

const $ = (id) => document.getElementById(id);

// Build the checkbox grid once. Subsequent scheme changes call
// refreshChainAvailability() to disable/enable checkboxes in place.
function buildChainCheckboxes() {
  const container = $("chain");
  container.innerHTML = "";
  for (const opt of ALL_CHAINS) {
    const label = document.createElement("label");
    label.className = "immunum-chain-check";
    const cb = document.createElement("input");
    cb.type = "checkbox";
    cb.value = opt.value;
    cb.checked = true;
    cb.dataset.chain = opt.value;
    const span = document.createElement("span");
    span.textContent = opt.label;
    label.appendChild(cb);
    label.appendChild(span);
    container.appendChild(label);
  }
}

function getScheme() {
  const el = document.querySelector('input[name="scheme"]:checked');
  return el ? el.value : "imgt";
}

function refreshChainAvailability() {
  const scheme = getScheme();
  const boxes = $("chain").querySelectorAll('input[type="checkbox"]');
  for (const cb of boxes) {
    const allowed = isChainAllowed(cb.value, scheme);
    cb.disabled = !allowed;
    cb.closest("label").classList.toggle("is-disabled", !allowed);
    if (!allowed) cb.checked = false;
  }
  // If the scheme change left everything unchecked, re-check all allowed
  // chains so the form stays submittable.
  const anyChecked = Array.from(boxes).some((cb) => cb.checked);
  if (!anyChecked) {
    for (const cb of boxes) {
      if (!cb.disabled) cb.checked = true;
    }
  }
}

function getSelectedChains() {
  return Array.from(
    $("chain").querySelectorAll('input[type="checkbox"]:checked'),
  ).map((cb) => cb.value);
}

function showError(msg) {
  const el = $("immunum-error");
  el.textContent = msg;
  el.hidden = false;
  $("immunum-result").hidden = true;
}

function clearError() {
  $("immunum-error").hidden = true;
}

// Map each index in the aligned sequence to a region label, based on the
// segment() substrings. Segments are contiguous and in order.
function buildRegionMap(sequence, segments) {
  const map = new Map();
  let cursor = 0;
  for (const region of REGIONS) {
    const seg = segments[region] || "";
    if (!seg.length) continue;
    // Find the segment starting at or after cursor (segments are non-overlapping
    // and appear in order in the original sequence).
    const idx = sequence.indexOf(seg, cursor);
    if (idx < 0) continue;
    for (let i = 0; i < seg.length; i++) {
      map.set(idx + i, region);
    }
    cursor = idx + seg.length;
  }
  return map;
}

function renderResult(sequence, numberResult, segResult) {
  const chainEl = $("result-chain");
  chainEl.textContent = `Chain: ${numberResult.chain}`;
  const schemeEl = $("result-scheme");
  schemeEl.textContent = `Scheme: ${numberResult.scheme}`;

  // Aligned query range (1-indexed, inclusive — user-facing convention).
  const qStart = numberResult.query_start;
  const qEnd = numberResult.query_end;
  $("result-range").textContent = `Query ${qStart + 1}–${qEnd + 1} (${qEnd - qStart + 1} aa)`;

  // Render the aligned region with the pre/post context dimmed.
  const alignedEl = $("result-aligned");
  alignedEl.textContent = "";
  const prefix = sequence.slice(0, qStart);
  const aligned = sequence.slice(qStart, qEnd + 1);
  const suffix = sequence.slice(qEnd + 1);
  if (prefix) {
    const s = document.createElement("span");
    s.className = "immunum-aligned-flank";
    s.textContent = prefix;
    alignedEl.appendChild(s);
  }
  const mid = document.createElement("span");
  mid.className = "immunum-aligned-core";
  mid.textContent = aligned;
  alignedEl.appendChild(mid);
  if (suffix) {
    const s = document.createElement("span");
    s.className = "immunum-aligned-flank";
    s.textContent = suffix;
    alignedEl.appendChild(s);
  }

  const conf = Math.max(0, Math.min(1, numberResult.confidence));
  $("result-confidence-value").textContent = `Confidence: ${(conf * 100).toFixed(1)}%`;

  // Figure out which region each numbered residue belongs to by locating the
  // aligned residue in the original sequence. We consume characters left to
  // right so repeated residues map in order.
  const regionMap = buildRegionMap(sequence, segResult);

  // Walk the sequence and the ordered numbering entries together. The WASM
  // bindings emit `numbering` by iterating positions in order, and JS object
  // iteration preserves insertion order for string keys.
  const entries = Array.from(numberResult.numbering);
  const grid = $("result-grid");
  grid.innerHTML = "";

  let seqCursor = 0;
  for (const [pos, aa] of entries) {
    const idx = sequence.indexOf(aa, seqCursor);
    let region = null;
    if (idx >= 0) {
      region = regionMap.get(idx) || null;
      seqCursor = idx + 1;
    }
    const pill = document.createElement("div");
    pill.className = "immunum-residue" + (region ? ` region-${region}` : "");
    pill.innerHTML = `<span class="pos">${pos}</span><span class="aa">${aa}</span>`;
    grid.appendChild(pill);
  }

  $("immunum-result").hidden = false;
}

async function main() {
  // Wire UI immediately so example buttons + scheme switching work even
  // before the WASM module finishes loading.
  buildChainCheckboxes();
  refreshChainAvailability();

  for (const r of document.querySelectorAll('input[name="scheme"]')) {
    r.addEventListener("change", refreshChainAvailability);
  }

  for (const btn of document.querySelectorAll("[data-example]")) {
    btn.addEventListener("click", () => {
      const key = btn.dataset.example;
      $("seq-input").value = EXAMPLES[key];
    });
  }

  try {
    await init();
  } catch (err) {
    showError(`Failed to load WASM module: ${err}`);
    return;
  }

  const runBtn = $("run-btn");
  runBtn.disabled = false;
  runBtn.textContent = "Number sequence";

  $("immunum-form").addEventListener("submit", (ev) => {
    ev.preventDefault();
    clearError();
    const sequence = $("seq-input").value.trim().toUpperCase().replace(/\s+/g, "");
    if (!sequence) {
      showError("Please enter a sequence.");
      return;
    }
    const scheme = getScheme();
    const chains = getSelectedChains();
    if (chains.length === 0) {
      showError("Please select at least one chain.");
      return;
    }

    let annotator;
    try {
      annotator = new Annotator(chains, scheme);
    } catch (err) {
      showError(`Could not initialise annotator: ${err}`);
      return;
    }
    try {
      const numResult = annotator.number(sequence);
      const segResult = annotator.segment(sequence);
      renderResult(sequence, numResult, segResult);
    } catch (err) {
      showError(`Numbering failed: ${err}`);
    } finally {
      annotator.free?.();
    }
  });
}

main();
