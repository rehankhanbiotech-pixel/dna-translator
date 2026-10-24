const test = require("node:test");
const assert = require("node:assert/strict");
const C = require("../dna-core.js");
const ok = (s) => { const r = C.parse(s); assert.ok(r.ok, r.error); return r.seq; };

test("A: ATGGCTTAA", () => {
  const s = ok("ATGGCTTAA");
  assert.equal(C.transcribe(s), "AUGGCUUAA");
  assert.equal(C.revComp(s), "TTAAGCCAT");
  assert.equal(C.translateFrame(s, "+1").protein, "MA*");
  assert.ok(Math.abs(C.stats(s).gc - 33.3333) < 0.01);
});
test("B: nested ORFs in +1", () => {
  const o = C.findOrfs("ATGATGTAA", 1);
  assert.equal(o.length, 2);
  assert.deepEqual(o.map((x) => [x.frame, x.protein, x.aa, x.from, x.to, x.complete]),
    [["+1", "MM*", 2, 1, 9, true], ["+1", "M*", 1, 4, 9, true]]);
});
test("C: CCCCCC has no ORF", () => assert.equal(C.findOrfs("CCCCCC", 1).length, 0));
test("D: ATGG ignores incomplete codon", () => {
  assert.equal(C.translateFrame("ATGG", "+1").protein, "M");
  const o = C.findOrfs("ATGG", 1);
  assert.equal(o.length, 1); assert.equal(o[0].complete, false); assert.equal(o[0].stopCodon, false);
});
test("E: N is rejected, never silently removed", () => {
  const r = C.parse("ATGNNNTAA");
  assert.equal(r.ok, false); assert.match(r.error, /Ambiguous base 'N' at position 4/);
});
test("F: RNA input equals DNA result", () => {
  const r = C.parse("AUGGCUUAA");
  assert.ok(r.ok && r.rna); assert.equal(r.seq, "ATGGCTTAA");
  assert.equal(C.translateFrame(r.seq, "+1").protein, "MA*");
  assert.equal(C.parse("ATGU").ok, false);
});
test("G: reverse-strand ORF with original coordinates", () => {
  const s = "CCTTAAGCCATG"; // contains rc of ATGGCTTAA at 3..11
  const o = C.findOrfs(s, 1).filter((x) => x.strand === "-");
  assert.equal(o.length, 1);
  assert.deepEqual([o[0].frame, o[0].protein, o[0].from, o[0].to, o[0].dna], ["-2", "MA*", 3, 11, "ATGGCTTAA"]);
  assert.equal(s.slice(2, 11), "TTAAGCCAT");
});
test("H: minimum ORF length counts amino acids before stop", () => {
  assert.deepEqual(C.findOrfs("ATGATGTAA", 2).map((x) => x.protein), ["MM*"]);
  assert.equal(C.findOrfs("ATGATGTAA", 3).length, 0);
});
test("input handling", () => {
  assert.match(C.parse("").error, /empty/i);
  assert.match(C.parse(">only header\n").error, /No nucleotides/);
  assert.match(C.parse("ATG1").error, /Unsupported character '1'/);
  assert.equal(ok("atggcttaa"), "ATGGCTTAA");
  const f = C.parse(">seq1 test\r\nATGG\r\n CT TAA\n\n");
  assert.ok(f.ok); assert.equal(f.seq, "ATGGCTTAA"); assert.equal(f.header, "seq1 test");
  assert.equal(C.parse(">a\nATG\n>b\nCCC").warnings.length, 1);
});
test("no stop codon -> incomplete ORF; stats never NaN", () => {
  const o = C.findOrfs("ATGGCTGCT", 1);
  assert.equal(o.length, 1); assert.equal(o[0].complete, false); assert.equal(o[0].protein, "MAA");
  const e = C.stats("");
  assert.equal(e.gc, 0); assert.ok(Number.isFinite(e.gc)); assert.equal(e.length, 0);
  assert.deepEqual(C.aaComposition("MA*A"), { counts: { M: 1, A: 2 }, length: 3 });
});
test("six frames", () => {
  assert.equal(C.translateFrame("AATGGCTTAA", "+2").protein, "MA*");
  assert.equal(C.translateFrame("ATGGCTTAA", "-1").protein, "LSH"); // TTA AGC CAT
  assert.equal(C.translateFrame("A", "+1").protein, "");
  assert.equal(Object.keys(C.CODON_TABLE).length, 64);
});
test("random sequences: coordinates round-trip and ORFs match a naive scan", () => {
  let seed = 7; const rnd = () => (seed = (seed * 1103515245 + 12345) % 2147483648) / 2147483648;
  for (let t = 0; t < 200; t++) {
    let s = ""; for (let i = 0, n = 5 + Math.floor(rnd() * 80); i < n; i++) s += "ACGT"[Math.floor(rnd() * 4)];
    const orfs = C.findOrfs(s, 1), L = s.length; let expected = 0;
    for (const strand of ["+", "-"]) {
      const t2 = strand === "+" ? s : C.revComp(s);
      for (let p = 0; p + 3 <= L; p++) if (t2.substr(p, 3) === "ATG") expected++;
    }
    assert.equal(orfs.length, expected); // one ORF per ATG, any frame
    for (const o of orfs) {
      const seg = s.slice(o.from - 1, o.to);
      assert.equal(o.strand === "+" ? seg : C.revComp(seg), o.dna);
      assert.equal(o.nt % 3, 0);
      if (o.complete) assert.equal(o.protein.at(-1), "*");
    }
  }
});

