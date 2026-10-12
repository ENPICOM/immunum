// Interactive WASM numbering tool for the immunum docs demo page.
import init, { Annotator, schemeSupportsChain, version } from "./wasm/immunum.js";

const EXAMPLES = {
  IGH: "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS",
  IGL: "SALTQPPAVSGTPGQRVTISCSGSDIGRRSVNWYQQFPGTAPKLLIYSNDQRPSVVPDRFSGSKSGTSASLAISGLQSEDEAEYYCAAWDDSLAVFGGGTQLTVGQPKA",
  TRB: "GVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASRPGLAGGRPEQYFGPGTRLTVTE",
  // The IGH example, a (G4S)3 linker and a kappa light chain.
  scFv:
    "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS" +
    "GGGGSGGGGSGGGGS" +
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK",
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

const $ = (id) => document.getElementById(id);

function el(tag, className, text) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
}

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

// Needs the WASM module: which chains a scheme numbers comes from immunum itself.
function refreshChainAvailability() {
  const scheme = getScheme();
  const boxes = $("chain").querySelectorAll('input[type="checkbox"]');
  for (const cb of boxes) {
    const allowed = schemeSupportsChain(scheme, cb.value);
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

// A segment result's regions as [name, residues] pairs, in sequence order and flanks included:
// segmentDomains() sets them in that order, and `error`/`errorKind` are its only other fields.
function segments(segResult) {
  return Object.entries(segResult).filter(([, residues]) => typeof residues === "string");
}

// "fr1" → "FR1", "prefix" → "Prefix".
function segmentLabel(name) {
  return /^(fr|cdr)\d$/.test(name) ? name.toUpperCase() : name[0].toUpperCase() + name.slice(1);
}

// The region of every residue in the sequence. All domains' segments in order rebuild the
// sequence, so this covers it end to end.
function residueRegions(segResults) {
  const regions = [];
  for (const segResult of segResults) {
    for (const [name, residues] of segments(segResult)) {
      for (let i = 0; i < residues.length; i++) regions.push(name);
    }
  }
  return regions;
}

// The whole sequence, coloured by region; flanks and linkers between domains are dimmed.
function renderSequence(segResults) {
  const pre = $("result-sequence");
  pre.textContent = "";
  segResults.forEach((segResult, d) => {
    for (const [name, residues] of segments(segResult)) {
      if (!residues) continue;
      const span = el("span", `region-${name}`, residues);
      span.title = `Domain ${d + 1}: ${segmentLabel(name)} (${residues.length} aa)`;
      pre.appendChild(span);
    }
  });
}

function renderSegmentTable(segResult) {
  const table = el("table", "immunum-segments");
  const head = el("tr");
  for (const title of ["Region", "Residues", "Length"]) head.appendChild(el("th", "", title));
  table.appendChild(el("thead")).appendChild(head);
  const body = table.appendChild(el("tbody"));
  for (const [name, residues] of segments(segResult)) {
    // Most domains have no flanks; an empty FR or CDR row still says something.
    if (!residues && (name === "prefix" || name === "postfix")) continue;
    const row = el("tr");
    const label = el("td");
    label.appendChild(el("span", `immunum-region-tag region-${name}`, segmentLabel(name)));
    row.append(label, el("td", "immunum-segment-residues", residues), el("td", "", residues.length));
    body.appendChild(row);
  }
  return table;
}

function renderDomain(index, numberResult, segResult, regions) {
  const section = el("section", "immunum-domain");

  // Aligned query range (1-indexed, inclusive — user-facing convention).
  const qStart = numberResult.queryStart;
  const qEnd = numberResult.queryEnd;
  const conf = Math.max(0, Math.min(1, numberResult.confidence));
  const header = el("div", "immunum-result-header");
  header.append(
    el("span", "immunum-domain-title", `Domain ${index + 1}`),
    el("span", "immunum-badge", `Chain: ${numberResult.chain}`),
    el("span", "immunum-badge", `Query ${qStart + 1}–${qEnd + 1} (${qEnd - qStart + 1} aa)`),
    el("span", "immunum-badge", `Confidence: ${(conf * 100).toFixed(1)}%`),
  );
  section.appendChild(header);

  section.appendChild(el("span", "immunum-aligned-label", "Segments"));
  section.appendChild(renderSegmentTable(segResult));

  section.appendChild(el("span", "immunum-aligned-label", "Numbering"));
  const grid = section.appendChild(el("div", "immunum-grid"));
  // The numbered residues run from queryStart to queryEnd, in whole-sequence coordinates.
  let seqIndex = qStart;
  for (const [pos, aa] of numberResult.numbering) {
    const region = regions[seqIndex++];
    const pill = el("div", "immunum-residue" + (region ? ` region-${region}` : ""));
    pill.append(el("span", "pos", pos), el("span", "aa", aa));
    grid.appendChild(pill);
  }
  return section;
}

// numberDomains() and segmentDomains() find the same domains, so their results pair up by index.
function renderResult(numberResults, segResults) {
  const count = numberResults.length;
  $("result-domain-count").textContent = `${count} domain${count === 1 ? "" : "s"}`;
  $("result-scheme").textContent = `Scheme: ${numberResults[0].scheme}`;
  renderSequence(segResults);

  const regions = residueRegions(segResults);
  const list = $("result-domains");
  list.textContent = "";
  numberResults.forEach((numberResult, d) => {
    list.appendChild(renderDomain(d, numberResult, segResults[d], regions));
  });

  $("immunum-result").hidden = false;
}

let wasmReady = null;

async function main() {
  // The script may run on a page that doesn't have the tool (e.g. when
  // Material's instant navigation loads a different page). Bail out early.
  if (!document.getElementById("immunum-form")) return;

  // Wire UI immediately so example buttons work even before the WASM
  // module finishes loading.
  buildChainCheckboxes();

  const confSlider = $("min-confidence");
  const confValue = $("confidence-value");
  confSlider.addEventListener("input", () => {
    confValue.textContent = parseFloat(confSlider.value).toFixed(2);
  });

  for (const btn of document.querySelectorAll("[data-example]")) {
    btn.addEventListener("click", () => {
      const key = btn.dataset.example;
      $("seq-input").value = EXAMPLES[key];
    });
  }

  try {
    if (!wasmReady) wasmReady = init();
    await wasmReady;
  } catch (err) {
    wasmReady = null;
    showError(`Failed to load WASM module: ${err}`);
    return;
  }

  // Read from the module itself, so the page can't claim a version it isn't running.
  const versionEl = $("immunum-version");
  versionEl.textContent = `Running immunum v${version()} (WebAssembly)`;
  versionEl.hidden = false;

  refreshChainAvailability();
  for (const r of document.querySelectorAll('input[name="scheme"]')) {
    r.addEventListener("change", refreshChainAvailability);
  }

  const runBtn = $("run-btn");
  runBtn.disabled = false;
  runBtn.textContent = "Number sequence";

  $("immunum-form").addEventListener("submit", (ev) => {
    ev.preventDefault();
    clearError();
    // A pasted sequence often wraps over lines or comes in spaced blocks; drop the whitespace, as
    // the CLI joins FASTA lines. Case is left alone: immunum numbers either case and returns
    // residues as given, like every interface.
    const sequence = $("seq-input").value.replace(/\s+/g, "");
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

    const minConfidence = parseFloat($("min-confidence").value);

    let annotator;
    try {
      annotator = new Annotator(chains, scheme, minConfidence);
    } catch (err) {
      showError(`Could not initialise annotator: ${err.message}`);
      return;
    }
    try {
      // Never empty: a sequence without a domain, or an invalid one, is a single result with
      // `error` set.
      const numResults = annotator.numberDomains(sequence);
      if (numResults[0].error) {
        showError(`Numbering failed: ${numResults[0].error}.`);
        return;
      }
      const segResults = annotator.segmentDomains(sequence);
      renderResult(numResults, segResults);
    } catch (err) {
      showError(`Numbering failed: ${err.message}`);
    } finally {
      annotator.free?.();
    }
  });
}

// Material for MkDocs uses instant navigation, which replaces the page
// body without re-running <script> tags. Subscribe to its `document$`
// observable so main() runs on every navigation to this page. Fall back
// to a plain invocation if the observable isn't available.
if (typeof window !== "undefined" && window.document$?.subscribe) {
  window.document$.subscribe(() => {
    main();
  });
} else {
  main();
}
