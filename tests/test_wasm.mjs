/**
 * WASM integration tests
 *
 * Prerequisites:
 *   wasm-pack build --target nodejs --features wasm --no-default-features
 *
 * Run:
 *   node --test tests/test_wasm.mjs
 */

import { strict as assert } from "node:assert";
import { readFileSync } from "node:fs";
import { describe, it } from "node:test";
import { Annotator, regionsFor, schemeSupportsChain } from "../pkg/immunum.js";

// The errors every interface reports, with their kind and message
const ERROR_CASES = JSON.parse(readFileSync(new URL("./error_cases.json", import.meta.url)));

// Asserts that `fn` throws an `Error` with the `kind` and `message` of `expected`
function assertThrowsCase(fn, expected) {
  assert.throws(fn, (err) => {
    assert.ok(err instanceof Error, `expected an Error, got ${err}`);
    assert.equal(err.kind, expected.kind);
    assert.equal(err.message, expected.message);
    return true;
  });
}

const ALL_CHAINS = ["H", "K", "L", "A", "B", "G", "D"];
const AB_CHAINS = ["H", "K", "L"];

const IGH_SEQ =
  "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";
const IGL_SEQ =
  "SALTQPPAVSGTPGQRVTISCSGSDIGRRSVNWYQQFPGTAPKLLIYSNDQRPSVVPDRFSGSKSGTSASLAISGLQSEDEAEYYCAAWDDSLAVFGGGTQLTVGQPKA";
const TRB_SEQ =
  "GVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASRPGLAGGRPEQYFGPGTRLTVTE";

describe("Annotator init", () => {
  it("constructs with short-form chains", () => {
    const annotator = new Annotator(["H"], "imgt");
    assert.ok(annotator);
  });

  it("constructs with long-form chains (IGH)", () => {
    const annotator = new Annotator(["IGH"], "imgt");
    assert.ok(annotator);
  });

  it("constructs with name-form chains (heavy)", () => {
    const annotator = new Annotator(["heavy"], "imgt");
    assert.ok(annotator);
  });

  for (const c of ERROR_CASES.setup) {
    it(`throws ${c.kind} for chains ${JSON.stringify(c.chains)}, scheme ${c.scheme}, minConfidence ${c.min_confidence}`, () => {
      assertThrowsCase(() => new Annotator(c.chains, c.scheme, c.min_confidence ?? undefined), c);
    });
  }

  it("throws unsupported_chain on antibody-only scheme + TCR", () => {
    for (const scheme of ["kabat", "chothia", "martin", "aho"]) {
      assert.throws(() => new Annotator(["A"], scheme), { kind: "unsupported_chain" });
    }
  });

  it("accepts every scheme name and its short alias", () => {
    for (const [scheme, alias, canonical] of [
      ["imgt", "i", "IMGT"],
      ["kabat", "k", "Kabat"],
      ["chothia", "c", "Chothia"],
      ["martin", "m", "Martin"],
      ["aho", "a", "Aho"],
    ]) {
      const byName = new Annotator(AB_CHAINS, scheme).number(IGH_SEQ);
      const byAlias = new Annotator(AB_CHAINS, alias).number(IGH_SEQ);
      assert.equal(byName.scheme, canonical);
      assert.equal(byAlias.scheme, canonical);
      assert.deepEqual([...byAlias.numbering], [...byName.numbering]);
    }
  });

  it("accepts chain groups like every other interface", () => {
    const byGroup = new Annotator(["ig"], "imgt").number(IGH_SEQ);
    const byChains = new Annotator(AB_CHAINS, "imgt").number(IGH_SEQ);
    assert.deepEqual([...byGroup.numbering], [...byChains.numbering]);
  });

  it("throws invalid_min_confidence below 0", () => {
    assert.throws(() => new Annotator(["H"], "imgt", -0.1), { kind: "invalid_min_confidence" });
  });
});

describe("schemeSupportsChain()", () => {
  it("allows every chain under IMGT and only antibody chains otherwise", () => {
    for (const chain of ALL_CHAINS) {
      assert.equal(schemeSupportsChain("imgt", chain), true);
    }
    for (const scheme of ["kabat", "chothia", "martin", "aho"]) {
      for (const chain of ALL_CHAINS) {
        assert.equal(schemeSupportsChain(scheme, chain), AB_CHAINS.includes(chain));
      }
    }
  });

  it("agrees with what the Annotator constructor accepts", () => {
    for (const scheme of ["imgt", "kabat", "chothia", "martin", "aho"]) {
      for (const chain of ALL_CHAINS) {
        let constructs = true;
        try {
          new Annotator([chain], scheme).free();
        } catch {
          constructs = false;
        }
        assert.equal(schemeSupportsChain(scheme, chain), constructs);
      }
    }
  });

  for (const c of ERROR_CASES.lookups) {
    const expected = c.scheme_supports_chain;
    it(`${typeof expected === "boolean" ? `returns ${expected}` : `throws ${expected.kind}`} for ${c.scheme}, ${c.chain}`, () => {
      if (typeof expected === "boolean") {
        assert.equal(schemeSupportsChain(c.scheme, c.chain), expected);
      } else {
        assertThrowsCase(() => schemeSupportsChain(c.scheme, c.chain), expected);
      }
    });
  }
});

describe("regionsFor()", () => {
  it("returns inclusive bounds in N- to C-terminal order", () => {
    const kabatHeavy = regionsFor("kabat", "H");
    assert.deepEqual(Object.keys(kabatHeavy), ["fr1", "cdr1", "fr2", "cdr2", "fr3", "cdr3", "fr4"]);
    assert.deepEqual(kabatHeavy.cdr1, [31, 35]);
    assert.deepEqual(kabatHeavy.fr4, [103, 113]);
    assert.notDeepEqual(regionsFor("kabat", "K"), kabatHeavy);
  });

  it("throws exactly for the pairs schemeSupportsChain rejects", () => {
    for (const scheme of ["imgt", "kabat", "chothia", "martin", "aho"]) {
      for (const chain of ALL_CHAINS) {
        if (schemeSupportsChain(scheme, chain)) {
          assert.doesNotThrow(() => regionsFor(scheme, chain));
        } else {
          assert.throws(() => regionsFor(scheme, chain), { kind: "unsupported_chain" });
        }
      }
    }
  });

  for (const c of ERROR_CASES.lookups) {
    it(`throws ${c.regions_for.kind} for ${c.scheme}, ${c.chain}`, () => {
      assertThrowsCase(() => regionsFor(c.scheme, c.chain), c.regions_for);
    });
  }
});

describe("number()", () => {
  it("returns correct chain and scheme for IGH", () => {
    const annotator = new Annotator(ALL_CHAINS, "imgt");
    const result = annotator.number(IGH_SEQ);
    assert.equal(result.chain, "H");
    assert.equal(result.scheme, "IMGT");
  });

  it("returns confidence as a number between 0 and 1", () => {
    const annotator = new Annotator(["H"], "imgt");
    const result = annotator.number(IGH_SEQ);
    assert.equal(typeof result.confidence, "number");
    assert.ok(result.confidence > 0 && result.confidence <= 1);
  });

  it("returns numbering as an ordered Map of string→residue", () => {
    const annotator = new Annotator(["H"], "imgt");
    const result = annotator.number(IGH_SEQ);
    assert.ok(result.numbering instanceof Map);
    assert.ok(result.numbering.size > 0);
    for (const [pos, aa] of result.numbering) {
      assert.equal(typeof pos, "string");
      assert.equal(typeof aa, "string");
      assert.equal(aa.length, 1);
    }
  });

  it("iterates numbering in IMGT order with insertions between base positions", () => {
    const annotator = new Annotator(["H"], "imgt");
    const result = annotator.number(IGH_SEQ);
    const keys = Array.from(result.numbering.keys());
    const idx111 = keys.indexOf("111");
    assert.ok(idx111 >= 0, "expected base position 111 to be present");
    const insertionKeys = keys.filter((k) => /^11[12][A-Z]$/.test(k));
    for (const ins of insertionKeys) {
      const insIdx = keys.indexOf(ins);
      assert.ok(insIdx > idx111, `insertion ${ins} must appear after 111`);
      const after = keys.slice(insIdx + 1);
      for (const k of after) {
        const m = /^(\d+)$/.exec(k);
        if (m) {
          assert.ok(
            parseInt(m[1], 10) >= 112,
            `insertion ${ins} should not precede base position ${k}`,
          );
          break;
        }
      }
    }
  });

  it("returns queryStart/queryEnd as inclusive 0-indexed ints on success", () => {
    const annotator = new Annotator(["H"], "imgt");
    const result = annotator.number(IGH_SEQ);
    assert.equal(typeof result.queryStart, "number");
    assert.equal(typeof result.queryEnd, "number");
    assert.ok(Number.isInteger(result.queryStart));
    assert.ok(Number.isInteger(result.queryEnd));
    assert.ok(result.queryStart >= 0);
    assert.ok(result.queryEnd >= result.queryStart);
    assert.ok(result.queryEnd < IGH_SEQ.length);
    // Aligned length should match the number of numbered residues.
    assert.equal(
      result.queryEnd - result.queryStart + 1,
      result.numbering.size,
    );
  });

  for (const c of ERROR_CASES.sequences) {
    it(`returns ${c.kind} with every other field null for ${c.sequence} (does not throw)`, () => {
      const result = new Annotator(c.chains, c.scheme).number(c.sequence);
      assert.deepEqual(result, {
        chain: null,
        scheme: null,
        confidence: null,
        numbering: null,
        queryStart: null,
        queryEnd: null,
        error: c.message,
        errorKind: c.kind,
      });
    });
  }

  it("returns null error and errorKind on success", () => {
    const annotator = new Annotator(["H"], "imgt");
    const result = annotator.number(IGH_SEQ);
    assert.equal(result.error, null);
    assert.equal(result.errorKind, null);
  });

  it("detects kappa or lambda for IGL sequence", () => {
    const annotator = new Annotator(AB_CHAINS, "imgt");
    const result = annotator.number(IGL_SEQ);
    assert.ok(["K", "L"].includes(result.chain));
  });

  it("returns kabat scheme label for kabat annotator", () => {
    const annotator = new Annotator(["H"], "kabat");
    const result = annotator.number(IGH_SEQ);
    assert.equal(result.scheme, "Kabat");
  });
});

describe("numberDomains() and segmentDomains()", () => {
  const LINKER = "GGGGSGGGGSGGGGS";
  // Starts at IMGT position 1. A light chain missing its first positions takes linker residues for
  // them when it follows a linker, so it wouldn't number as it does on its own.
  const KAPPA_SEQ =
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK";
  const REGIONS = ["prefix", "fr1", "cdr1", "fr2", "cdr2", "fr3", "cdr3", "fr4", "postfix"];

  it("numbers each domain like number() on its own", () => {
    const annotator = new Annotator(["ig"], "imgt");
    const domains = annotator.numberDomains(IGH_SEQ + LINKER + KAPPA_SEQ);
    const heavy = annotator.number(IGH_SEQ);
    const light = annotator.number(KAPPA_SEQ);
    assert.deepEqual(
      domains.map((d) => d.chain),
      [heavy.chain, light.chain],
    );
    assert.deepEqual([...domains[0].numbering], [...heavy.numbering]);
    assert.deepEqual([...domains[1].numbering], [...light.numbering]);
    assert.equal(domains[1].queryStart, IGH_SEQ.length + LINKER.length + light.queryStart);
  });

  it("segments with every residue in exactly one domain", () => {
    const annotator = new Annotator(["ig"], "imgt");
    const sequence = "MKYLL" + IGH_SEQ + LINKER + KAPPA_SEQ + "HHHHHH";
    const domains = annotator.segmentDomains(sequence);
    assert.equal(domains.map((d) => REGIONS.map((r) => d[r]).join("")).join(""), sequence);
    assert.deepEqual(
      domains.map((d) => [d.prefix, d.postfix]),
      [
        ["MKYLL", ""],
        [LINKER, "HHHHHH"],
      ],
    );
    assert.equal(domains[1].cdr3, annotator.segment(KAPPA_SEQ).cdr3);
  });

  for (const method of ["numberDomains", "segmentDomains"]) {
    for (const c of [...ERROR_CASES.sequences, ...ERROR_CASES.domain_errors]) {
      it(`${method} returns one ${c.kind} result for ${c.sequence}`, () => {
        const domains = new Annotator(c.chains, c.scheme)[method](c.sequence);
        assert.equal(domains.length, 1);
        assert.equal(domains[0].error, c.message);
        assert.equal(domains[0].errorKind, c.kind);
      });
    }
  }

  for (const c of ERROR_CASES.domain_errors) {
    it(`number() and segment() succeed for ${c.sequence}, which holds no domain`, () => {
      const annotator = new Annotator(c.chains, c.scheme);
      assert.equal(annotator.number(c.sequence).error, null);
      assert.equal(annotator.segment(c.sequence).error, null);
    });
  }
});

describe("segment()", () => {
  it("returns all expected region keys", () => {
    const annotator = new Annotator(ALL_CHAINS, "IMGT");
    const result = annotator.segment(IGH_SEQ);
    const expected = ["fr1", "cdr1", "fr2", "cdr2", "fr3", "cdr3", "fr4", "prefix", "postfix"];
    for (const key of expected) {
      assert.ok(key in result, `missing key: ${key}`);
      assert.equal(typeof result[key], "string");
    }
  });

  it("fr1 is non-empty for a full sequence", () => {
    const annotator = new Annotator(["IGH"], "IMGT");
    const result = annotator.segment(IGH_SEQ);
    assert.ok(result.fr1.length > 0);
  });

  it("works for TCR sequences", () => {
    const annotator = new Annotator(["TRB"], "IMGT");
    const result = annotator.segment(TRB_SEQ);
    assert.ok(result.cdr3.length > 0);
  });

  it("returns null error and errorKind on success", () => {
    const annotator = new Annotator(["H"], "IMGT");
    const result = annotator.segment(IGH_SEQ);
    assert.equal(result.error, null);
    assert.equal(result.errorKind, null);
  });

  for (const c of ERROR_CASES.sequences) {
    it(`returns ${c.kind} with every region null for ${c.sequence} (does not throw)`, () => {
      const result = new Annotator(c.chains, c.scheme).segment(c.sequence);
      assert.equal(result.error, c.message);
      assert.equal(result.errorKind, c.kind);
      for (const key of ["prefix", "fr1", "cdr1", "fr2", "cdr2", "fr3", "cdr3", "fr4", "postfix"]) {
        assert.equal(result[key], null, `${key} should be null`);
      }
    });
  }

  it("keeps flanking residues in prefix and postfix (#58)", () => {
    const annotator = new Annotator(["H"], "IMGT");
    const sequence = "MGWSCIILFLVATATGVHSX" + IGH_SEQ + "HHHHHHEPEA";
    const result = annotator.segment(sequence);
    assert.equal(result.prefix, "MGWSCIILFLVATATGVHSX");
    assert.equal(result.postfix, "HHHHHHEPEA");
    const regions = ["prefix", "fr1", "cdr1", "fr2", "cdr2", "fr3", "cdr3", "fr4", "postfix"];
    assert.equal(regions.map((r) => result[r]).join(""), sequence);
  });
});
