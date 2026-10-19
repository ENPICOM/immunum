// Interactive WASM numbering tool for the immunum docs demo page.
import init, { Annotator, regionsFor, schemeSupportsChain, version } from "./wasm/immunum.js";

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

// The scheme radio values, in the order the page lists them, with their display names.
const SCHEMES = [
  { value: "imgt", label: "IMGT" },
  { value: "kabat", label: "Kabat" },
  { value: "chothia", label: "Chothia" },
  { value: "martin", label: "Martin" },
  { value: "aho", label: "AHo" },
];

// The FR/CDR regions N- to C-terminal, and the segments of a domain in sequence order: the regions
// between the flanks, as every interface names them.
const REGIONS = ["fr1", "cdr1", "fr2", "cdr2", "fr3", "cdr3", "fr4"];
const SEGMENTS = ["prefix", ...REGIONS, "postfix"];
const CDRS = ["cdr1", "cdr2", "cdr3"];

const DEFAULT_MIN_CONFIDENCE = "0.4";

const $ = (id) => document.getElementById(id);

function el(tag, className, text) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
}

// Build the checkbox grid once per page. Scheme changes call refreshChainAvailability() to
// disable/enable checkboxes in place.
function buildChainCheckboxes() {
  const container = $("chain");
  container.textContent = "";
  for (const opt of ALL_CHAINS) {
    const label = el("label", "immunum-chain-check");
    const cb = el("input");
    cb.type = "checkbox";
    cb.value = opt.value;
    cb.checked = true;
    label.append(cb, el("span", "", opt.label));
    container.appendChild(label);
  }
}

function chainBoxes() {
  return Array.from($("chain").querySelectorAll('input[type="checkbox"]'));
}

function getScheme() {
  const checked = document.querySelector('input[name="scheme"]:checked');
  return checked ? checked.value : "imgt";
}

// Needs the WASM module: which chains a scheme numbers comes from immunum itself. A chain the
// scheme doesn't number is disabled but keeps its checked state, so switching back to a scheme
// that numbers it restores the selection.
function refreshChainAvailability() {
  const scheme = getScheme();
  const boxes = chainBoxes();
  for (const cb of boxes) {
    const allowed = schemeSupportsChain(scheme, cb.value);
    cb.disabled = !allowed;
    cb.closest("label").classList.toggle("is-disabled", !allowed);
  }
  // If the scheme leaves no usable chain checked, check every usable one so the form stays
  // submittable.
  if (!boxes.some((cb) => cb.checked && !cb.disabled)) {
    for (const cb of boxes) if (!cb.disabled) cb.checked = true;
  }
}

function getSelectedChains() {
  return chainBoxes()
    .filter((cb) => cb.checked && !cb.disabled)
    .map((cb) => cb.value);
}

function getMinConfidence() {
  return parseFloat($("min-confidence").value);
}

function showConfidence() {
  $("confidence-value").textContent = getMinConfidence().toFixed(2);
}

// The error and notice regions are live regions that stay in the page and hide when empty, so
// screen readers announce each new message.
function showError(msg) {
  $("immunum-error").textContent = msg;
  $("immunum-notice").textContent = "";
  $("immunum-result").hidden = true;
  current = null;
}

function clearMessages() {
  $("immunum-error").textContent = "";
  $("immunum-notice").textContent = "";
  $("immunum-action-status").textContent = "";
}

// The records of a pasted sequence or FASTA text. Lines starting with ">" are FASTA headers that
// start a record named by the header's first word, as the CLI names it. Every other line adds to
// the current record with its whitespace dropped, so a sequence wrapped over lines or written in
// spaced blocks is joined; text before the first header is a record named "seq_1", as the CLI names
// a bare sequence. Case is left alone: immunum numbers either case and returns residues as given.
// Records without residues are dropped, as the CLI drops them.
function parseInput(text) {
  const records = [];
  let record = null;
  for (const line of text.split(/\r?\n/)) {
    const trimmed = line.trim();
    if (trimmed.startsWith(">")) {
      record = { id: trimmed.slice(1).trim().split(/\s+/)[0] || "unknown", sequence: "" };
      records.push(record);
    } else if (trimmed) {
      if (!record) {
        record = { id: "seq_1", sequence: "" };
        records.push(record);
      }
      record.sequence += trimmed.replace(/\s+/g, "");
    }
  }
  return records.filter((r) => r.sequence);
}

// ---------- Rendering ----------

// The segments of a domain as [name, residues] pairs in sequence order, flanks included.
function segments(segResult) {
  return SEGMENTS.map((name) => [name, segResult[name] ?? ""]);
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

function table(className, headings) {
  const node = el("table", className);
  const head = el("tr");
  for (const title of headings) head.appendChild(el("th", "", title));
  node.appendChild(el("thead")).appendChild(head);
  return [node, node.appendChild(el("tbody"))];
}

function regionTag(name) {
  return el("span", `immunum-region-tag region-${name}`, segmentLabel(name));
}

function renderSegmentTable(segResult) {
  const [node, body] = table("immunum-segments", ["Region", "Residues", "Length"]);
  for (const [name, residues] of segments(segResult)) {
    // Most domains have no flanks; an empty FR or CDR row still says something.
    if (!residues && (name === "prefix" || name === "postfix")) continue;
    const row = el("tr");
    const label = el("td");
    label.appendChild(regionTag(name));
    row.append(label, el("td", "immunum-segment-residues", residues), el("td", "", residues.length));
    body.appendChild(row);
  }
  return node;
}

// Every scheme's numbering and segments of the sequence, keyed by scheme value, each from its own
// Annotator over the selected chains that scheme numbers. Null for a scheme that numbers none of
// them. Domain detection doesn't depend on the scheme, but the chains a scheme allows can change
// which domains are found, so domains are matched up by chain and range, not index.
function numberUnderEveryScheme(sequence, chains, minConfidence) {
  const results = {};
  for (const { value: scheme } of SCHEMES) {
    const allowed = chains.filter((c) => schemeSupportsChain(scheme, c));
    if (allowed.length === 0) {
      results[scheme] = null;
      continue;
    }
    const annotator = new Annotator(allowed, scheme, minConfidence);
    try {
      const numbers = annotator.numberDomains(sequence);
      results[scheme] = numbers[0].error ? null : { numbers, segments: annotator.segmentDomains(sequence) };
    } finally {
      annotator.free();
    }
  }
  return results;
}

// The segments `schemeResult` gives the domain `numberResult` found: those of its domain of the
// same chain whose range overlaps it.
function matchingSegments(schemeResult, numberResult) {
  if (!schemeResult) return null;
  const index = schemeResult.numbers.findIndex(
    (n) =>
      n.chain === numberResult.chain &&
      n.queryStart <= numberResult.queryEnd &&
      n.queryEnd >= numberResult.queryStart,
  );
  return index === -1 ? null : schemeResult.segments[index];
}

// FR/CDR boundaries and CDR residues of one domain under every scheme that numbers its chain.
function renderComparison(numberResult, everyScheme, scheme) {
  const details = el("details", "immunum-compare");
  details.appendChild(el("summary", "", "Compare schemes"));
  const schemes = SCHEMES.filter((s) => schemeSupportsChain(s.value, numberResult.chain));

  details.appendChild(el("span", "immunum-aligned-label", `Region boundaries for chain ${numberResult.chain}`));
  const [bounds, boundsBody] = table("immunum-segments immunum-compare-table", [
    "Scheme",
    ...REGIONS.map(segmentLabel),
  ]);
  for (const s of schemes) {
    const spans = regionsFor(s.value, numberResult.chain);
    const row = el("tr", s.value === scheme ? "is-current" : "");
    row.appendChild(el("th", "", s.label)).scope = "row";
    for (const name of REGIONS) row.appendChild(el("td", "", `${spans[name][0]}–${spans[name][1]}`));
    boundsBody.appendChild(row);
  }
  details.appendChild(scrollBox(bounds));

  details.appendChild(el("span", "immunum-aligned-label", "CDR residues"));
  const [cdrs, cdrsBody] = table("immunum-segments immunum-compare-table", ["Scheme", ...CDRS.map(segmentLabel)]);
  for (const s of schemes) {
    const seg = matchingSegments(everyScheme[s.value], numberResult);
    const row = el("tr", s.value === scheme ? "is-current" : "");
    row.appendChild(el("th", "", s.label)).scope = "row";
    for (const name of CDRS) {
      const cell = el("td", "immunum-segment-residues", seg ? seg[name] : "not found");
      if (seg) cell.title = `${seg[name].length} aa`;
      row.appendChild(cell);
    }
    cdrsBody.appendChild(row);
  }
  details.appendChild(scrollBox(cdrs));
  return details;
}

// Wide tables scroll inside their own box instead of widening the page.
function scrollBox(child) {
  const box = el("div", "immunum-scroll");
  box.appendChild(child);
  return box;
}

function renderDomain(index, numberResult, segResult, regions, everyScheme, scheme) {
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
  section.appendChild(renderComparison(numberResult, everyScheme, scheme));

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
function renderResult({ numbers, segments: segResults, everyScheme, scheme }) {
  const count = numbers.length;
  $("result-domain-count").textContent = `${count} domain${count === 1 ? "" : "s"}`;
  $("result-scheme").textContent = `Scheme: ${numbers[0].scheme}`;
  renderSequence(segResults);

  const regions = residueRegions(segResults);
  const list = $("result-domains");
  list.textContent = "";
  numbers.forEach((numberResult, d) => {
    list.appendChild(renderDomain(d, numberResult, segResults[d], regions, everyScheme, scheme));
  });

  $("immunum-result").hidden = false;
}

// ---------- Export ----------

// TSV in the columns of the CLI's `immunum number --all-domains`, one row per numbered residue:
//   sequence_id, domain (0-based), chain, scheme, confidence (4 decimals), position, residue,
//   error, error_kind
// The demo only exports numbered domains, so error and error_kind are always empty.
function toTsv({ id, numbers }) {
  const lines = ["sequence_id\tdomain\tchain\tscheme\tconfidence\tposition\tresidue\terror\terror_kind"];
  numbers.forEach((n, d) => {
    const confidence = n.confidence.toFixed(4);
    for (const [pos, aa] of n.numbering) {
      lines.push(`${id}\t${d}\t${n.chain}\t${n.scheme}\t${confidence}\t${pos}\t${aa}\t\t`);
    }
  });
  return lines.join("\n") + "\n";
}

// JSON: an array with one object per domain, keyed as a record of the CLI's
// `immunum number --all-domains --format json` (sequence_id, domain, chain, scheme, confidence,
// numbering, query_start, query_end, error, error_kind), plus `segments`: the domain's prefix,
// fr1 … fr4 and postfix, as `immunum segment` names them.
function toJson({ id, numbers, segments: segResults }) {
  const records = numbers.map((n, d) => ({
    sequence_id: id,
    domain: d,
    chain: n.chain,
    scheme: n.scheme,
    confidence: n.confidence,
    numbering: n.numbering,
    query_start: n.queryStart,
    query_end: n.queryEnd,
    error: null,
    error_kind: null,
    segments: Object.fromEntries(segments(segResults[d])),
  }));
  return stringify(records) + "\n";
}

// JSON.stringify(value, null, 2), except that a Map is written as an object in its own order:
// as a plain object, the numbering would list integer-like positions ("1", "112") before
// insertions ("111A").
function stringify(value, indent = "") {
  const inner = indent + "  ";
  if (Array.isArray(value)) {
    if (value.length === 0) return "[]";
    return `[\n${value.map((v) => inner + stringify(v, inner)).join(",\n")}\n${indent}]`;
  }
  if (value instanceof Map || (value !== null && typeof value === "object")) {
    const entries = value instanceof Map ? [...value] : Object.entries(value);
    if (entries.length === 0) return "{}";
    const lines = entries.map(([k, v]) => `${inner}${JSON.stringify(k)}: ${stringify(v, inner)}`);
    return `{\n${lines.join(",\n")}\n${indent}}`;
  }
  return JSON.stringify(value);
}

function download(text, type, extension) {
  const url = URL.createObjectURL(new Blob([text], { type }));
  const link = el("a");
  link.href = url;
  link.download = `immunum-${current.id.replace(/[^\w.-]+/g, "_")}.${extension}`;
  document.body.appendChild(link);
  link.click();
  link.remove();
  // Revoke once the click has been handled.
  setTimeout(() => URL.revokeObjectURL(url), 0);
}

// `failed` says where else to get the text when the browser blocks the clipboard.
async function copy(text, done, failed) {
  const status = $("immunum-action-status");
  try {
    await navigator.clipboard.writeText(text);
    status.textContent = done;
  } catch {
    status.textContent = `Copying failed: the browser blocked clipboard access. ${failed}`;
  }
}

// ---------- Shareable link ----------

// The query string holds the sequence and settings: ?seq=…&scheme=kabat&chains=HKL&min_confidence=0.6
function shareUrl(sequence) {
  const url = new URL(window.location.href);
  url.hash = "";
  url.search = new URLSearchParams({
    seq: sequence,
    scheme: getScheme(),
    chains: getSelectedChains().join(""),
    min_confidence: String(getMinConfidence()),
  }).toString();
  return url.toString();
}

// Prefill the form from the query string. True when it holds a sequence to number.
function applyUrl() {
  const params = new URLSearchParams(window.location.search);
  const scheme = params.get("scheme");
  if (scheme && SCHEMES.some((s) => s.value === scheme)) {
    document.querySelector(`input[name="scheme"][value="${scheme}"]`).checked = true;
  }
  const chains = params.get("chains");
  if (chains) {
    for (const cb of chainBoxes()) cb.checked = chains.toUpperCase().includes(cb.value);
  }
  const minConfidence = parseFloat(params.get("min_confidence"));
  if (minConfidence >= 0 && minConfidence <= 1) $("min-confidence").value = String(minConfidence);
  showConfidence();
  const seq = params.get("seq");
  if (seq) $("seq-input").value = seq;
  return Boolean(seq);
}

// ---------- Running ----------

// The shown result: the numbered record, every domain's numbering and segments, the same under
// every scheme for the comparison, and the scheme used. Null when no result is shown.
let current = null;
// The WASM module's load, started by the first main(), and whether it has finished.
let wasmReady = null;
let wasmLoaded = false;

function run() {
  if (!wasmLoaded) return;
  clearMessages();
  const records = parseInput($("seq-input").value);
  if (records.length === 0) {
    showError("Please enter a sequence.");
    return;
  }
  const { id, sequence } = records[0];
  const scheme = getScheme();
  const chains = getSelectedChains();
  if (chains.length === 0) {
    showError("Please select at least one chain.");
    return;
  }
  const minConfidence = getMinConfidence();

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
    const numbers = annotator.numberDomains(sequence);
    if (numbers[0].error) {
      showError(`Numbering failed: ${numbers[0].error}.`);
      return;
    }
    const segResults = annotator.segmentDomains(sequence);
    const everyScheme = numberUnderEveryScheme(sequence, chains, minConfidence);
    current = { id, sequence, numbers, segments: segResults, everyScheme, scheme };
    renderResult(current);
  } catch (err) {
    showError(`Numbering failed: ${err.message}`);
    return;
  } finally {
    annotator.free();
  }

  if (records.length > 1) {
    $("immunum-notice").textContent =
      `Numbered only the first of ${records.length} FASTA records (${id}).`;
  }
  // Keep the address bar in step with the result, so it can be bookmarked or shared as is.
  history.replaceState(history.state, "", shareUrl(sequence));
}

// After a settings change, a shown result is out of date: number again with the new settings.
// That clears the result and shows the error if the input no longer numbers.
function rerunIfShown() {
  if (current) run();
}

// Event listeners go on the page's own elements, which Material's instant navigation replaces
// on every visit, so wiring once per form never doubles them.
function wire() {
  const form = $("immunum-form");
  if (form.dataset.wired) return false;
  form.dataset.wired = "true";

  buildChainCheckboxes();

  $("min-confidence").addEventListener("input", showConfidence);
  $("min-confidence").addEventListener("change", rerunIfShown);
  $("chain").addEventListener("change", rerunIfShown);
  for (const radio of document.querySelectorAll('input[name="scheme"]')) {
    radio.addEventListener("change", () => {
      refreshChainAvailability();
      rerunIfShown();
    });
  }

  for (const btn of document.querySelectorAll("[data-example]")) {
    btn.addEventListener("click", () => {
      $("seq-input").value = EXAMPLES[btn.dataset.example];
    });
  }

  form.addEventListener("submit", (ev) => {
    ev.preventDefault();
    run();
  });
  $("seq-input").addEventListener("keydown", (ev) => {
    if (ev.key === "Enter" && (ev.ctrlKey || ev.metaKey)) {
      ev.preventDefault();
      run();
    }
  });

  $("copy-tsv").addEventListener("click", () =>
    copy(toTsv(current), "Copied the numbering as TSV.", "Use Download TSV instead."),
  );
  $("download-tsv").addEventListener("click", () =>
    download(toTsv(current), "text/tab-separated-values", "tsv"),
  );
  $("download-json").addEventListener("click", () =>
    download(toJson(current), "application/json", "json"),
  );
  $("copy-link").addEventListener("click", () =>
    copy(shareUrl(current.sequence), "Copied a link to this result.", "The address bar holds the same link."),
  );
  return true;
}

async function main() {
  // The script may run on a page that doesn't have the tool (e.g. when
  // Material's instant navigation loads a different page). Bail out early.
  if (!$("immunum-form")) return;

  // Wire the UI before the WASM module loads, so example buttons work meanwhile. A form that is
  // already wired has had its setup run.
  if (!wire()) return;
  current = null;
  const runFromUrl = applyUrl();

  try {
    // One shared load: main() runs again on every visit, possibly before an earlier load ends.
    wasmReady ??= init();
    await wasmReady;
  } catch (err) {
    wasmReady = null;
    showError(`Failed to load WASM module: ${err}`);
    return;
  }
  wasmLoaded = true;

  // Read from the module itself, so the page can't claim a version it isn't running.
  const versionEl = $("immunum-version");
  versionEl.textContent = `Running immunum v${version()} (WebAssembly)`;
  versionEl.hidden = false;

  refreshChainAvailability();

  const runBtn = $("run-btn");
  runBtn.disabled = false;
  runBtn.textContent = "Number sequence";

  if (runFromUrl) run();
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
