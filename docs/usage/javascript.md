# JavaScript

The npm package is the immunum core compiled to WebAssembly, built for Node.js, with TypeScript
types. Install with `npm install immunum` ([Getting started](../getting-started.md#javascript));
every type and function is in the [API reference](../api/javascript.md).

Results use camelCase: `queryStart`, `queryEnd` and `errorKind` are `query_start`, `query_end` and
`error_kind` in the other interfaces. Every other field name is the same.

## Numbering and segmentation

```js
const { Annotator } = require("immunum");

// Chains, scheme and an optional minimum confidence, as in Python
const annotator = new Annotator(["ig"], "imgt");

const sequence =
  "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";

const result = annotator.number(sequence);
console.log(result.chain);                       // "H"
console.log(result.numbering.get("111A"));       // "K": a Map, in sequence order
console.log(result.queryStart, result.queryEnd); // 0 121

const segments = annotator.segment(sequence);
console.log(segments.cdr3); // "AREGTTGKPIGAFAH"

annotator.free();
```

An `Annotator` holds WebAssembly memory. The garbage collector frees it eventually, but not on any
schedule, so free it yourself when you're done: call `free()`, or declare it with `using` (explicit
resource management: Node.js 24 or newer, or TypeScript 5.2 or newer) to free it at the end of the
block:

```js
{
  using kabat = new Annotator(["H"], "kabat");
  kabat.segment(sequence).cdr1; // "TYALS"
}
```

Chain and scheme names, and what the minimum confidence does, are in
[Concepts](../concepts.md).

## Multiple domains

`numberDomains` and `segmentDomains` return one result per variable domain, such as both domains of
an scFv, in sequence order. Each is what `number` or `segment` returns for that domain, and the
residues between two domains are the second domain's `prefix`:

```js
const heavy =
  "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";
const kappa =
  "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK";
const scfv = heavy + "GGGGSGGGGSGGGGS" + kappa;

const annotator = new Annotator(["ig"], "imgt");
annotator.numberDomains(scfv).map((d) => d.chain);    // ["H", "K"]
annotator.segmentDomains(scfv).map((s) => s.prefix);  // ["", "GGGGSGGGGSGGGGS"]
annotator.free();
```

## Schemes and regions

`schemeSupportsChain(scheme, chain)` tells whether a scheme numbers a chain before you build an
annotator, e.g. to offer only valid choices in a UI. `regionsFor(scheme, chain)` returns the
scheme's FR/CDR boundaries for the chain, both ends inclusive, without numbering a sequence.
`version()` returns the version of immunum the module was built from:

```js
const { regionsFor, schemeSupportsChain, version } = require("immunum");

schemeSupportsChain("kabat", "H"); // true
schemeSupportsChain("kabat", "B"); // false: only IMGT numbers TCR chains
regionsFor("kabat", "H").cdr1;     // [31, 35]
version();                         // "2.0.0"
```

## Errors

A setup mistake throws an `Error` with a `kind`; a sequence that can't be numbered comes back with
`error` and `errorKind` set and every other field `null`. See [Errors](errors.md#javascript).
