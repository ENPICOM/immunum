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
import { describe, it } from "node:test";
import { Annotator } from "../pkg/immunum.js";

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

  it("throws on invalid chain", () => {
    assert.throws(() => new Annotator(["INVALID"], "imgt"));
  });

  it("throws on invalid scheme", () => {
    assert.throws(() => new Annotator(["H"], "INVALID"));
  });

  it("throws on Kabat + TCR", () => {
    assert.throws(() => new Annotator(["A"], "kabat"));
  });
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

  it("returns numbering as an object with string keys and single-char values", () => {
    const annotator = new Annotator(["H"], "imgt");
    const result = annotator.number(IGH_SEQ);
    assert.equal(typeof result.numbering, "object");
    const entries = Object.entries(result.numbering);
    assert.ok(entries.length > 0);
    for (const [pos, aa] of entries) {
      assert.equal(typeof pos, "string");
      assert.equal(typeof aa, "string");
      assert.equal(aa.length, 1);
    }
  });

  it("returns error field on empty sequence (does not throw)", () => {
    const annotator = new Annotator(ALL_CHAINS, "imgt");
    const result = annotator.number("");
    assert.equal(result.chain, null);
    assert.equal(result.numbering, null);
    assert.equal(typeof result.error, "string");
  });

  it("returns error field on invalid sequence (does not throw)", () => {
    const annotator = new Annotator(ALL_CHAINS, "imgt");
    const result = annotator.number("AAAAAAAAAAAAAAAA");
    assert.equal(result.chain, null);
    assert.equal(result.numbering, null);
    assert.equal(typeof result.error, "string");
  });

  it("returns null error on success", () => {
    const annotator = new Annotator(["H"], "imgt");
    const result = annotator.number(IGH_SEQ);
    assert.equal(result.error, null);
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

  it("returns null error on success", () => {
    const annotator = new Annotator(["H"], "IMGT");
    const result = annotator.segment(IGH_SEQ);
    assert.equal(result.error, null);
  });

  it("returns error field on invalid sequence (does not throw)", () => {
    const annotator = new Annotator(ALL_CHAINS, "IMGT");
    const result = annotator.segment("AAAAAAAAAAAAAAAA");
    assert.equal(typeof result.error, "string");
    assert.equal(result.fr1, undefined);
  });
});
