/* dna-core.js - pure DNA logic (no DOM). Loaded by index.html in the browser and by Node tests.
 * Conventions: input is the CODING (sense) strand, 5'->3'. Coordinates are 1-based, inclusive,
 * always on the ORIGINAL input sequence (reverse-strand ORFs: from < to, strand "-").
 * Frame -N = reverse-complement sequence read from its (N)th base (offset N-1). Standard genetic code only. */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.DNACore = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";
  var B = "TCAG", AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", CT = {};
  for (var i = 0; i < 4; i++) for (var j = 0; j < 4; j++) for (var k = 0; k < 4; k++) CT[B[i] + B[j] + B[k]] = AAS[16 * i + 4 * j + k];
  var FRAMES = ["+1", "+2", "+3", "-1", "-2", "-3"], COMP = { A: "T", T: "A", G: "C", C: "G" };

  function revComp(s) { var o = ""; for (var n = s.length - 1; n >= 0; n--) o += COMP[s[n]]; return o; }
  function transcribe(s) { return s.replace(/T/g, "U"); } // coding strand -> mRNA (same 5'->3' sequence)

  // Accepts DNA/RNA, FASTA (first record used), multiline, any case. Rejects anything but A,C,G,T,U (no silent deletion).
  function parse(raw) {
    var r = { ok: false, seq: "", header: "", rna: false, warnings: [], error: "" };
    if (typeof raw !== "string" || !raw.trim()) { r.error = "Input is empty. Paste a DNA sequence or FASTA."; return r; }
    var recs = [], cur = null;
    raw.split(/\r\n|\r|\n/).forEach(function (l) {
      var t = l.trim(); if (!t) return;
      if (t[0] === ">") { cur = { h: t.slice(1).trim(), b: [] }; recs.push(cur); }
      else { if (!cur) { cur = { h: "", b: [] }; recs.push(cur); } cur.b.push(t); }
    });
    if (!recs.length || !recs[0].b.length) { r.error = "No nucleotides found (only a FASTA header or blank lines)."; return r; }
    if (recs.length > 1) r.warnings.push("Multiple FASTA records found; only the first (" + (recs[0].h || "unnamed") + ") was used.");
    r.header = recs[0].h;
    var s = recs[0].b.join("").replace(/\s+/g, "").toUpperCase();
    var badAll = s.match(/[^ACGTU]/g);
    if (badAll) {
      var ch = badAll[0], pos = s.search(/[^ACGTU]/) + 1;
      r.error = (/[NRYKMSWBDHV]/.test(ch) ? "Ambiguous base '" : "Unsupported character '") + ch + "' at position " + pos +
        " (" + badAll.length + " unsupported in total). Only A, C, G, T (or U) are accepted; ambiguity codes such as N are not supported, so remove or replace them.";
      return r;
    }
    if (/U/.test(s) && /T/.test(s)) { r.error = "Sequence contains both T and U. Enter DNA (T) or RNA (U), not a mixture."; return r; }
    if (/U/.test(s)) { s = s.replace(/U/g, "T"); r.rna = true; r.warnings.push("RNA input detected: U was converted to T (the sequence is treated as the coding strand)."); }
    r.seq = s; r.ok = true; return r;
  }

  // Translate one frame; trailing incomplete codon (1-2 nt) is ignored. pos = 0-based offset in the frame's strand.
  function translateFrame(seq, frame) {
    var s = frame[0] === "-" ? revComp(seq) : seq, o = +frame[1] - 1, cs = [];
    for (var n = o; n + 3 <= s.length; n += 3) { var c = s.substr(n, 3); cs.push({ c: c, a: CT[c], pos: n }); }
    return { codons: cs, protein: cs.map(function (x) { return x.a; }).join("") };
  }

  // ORF = ATG ... first in-frame stop (stop included in DNA/coordinates, excluded from aa length).
  // Every ATG is scanned independently, so nested/overlapping ORFs are all reported.
  // Without a stop the ORF runs to the last complete codon and complete=false. minAA counts amino acids before the stop.
  function findOrfs(seq, minAA) {
    minAA = minAA > 0 ? Math.floor(minAA) : 1;
    var L = seq.length, res = [];
    FRAMES.forEach(function (f) {
      var cs = translateFrame(seq, f).codons, o = +f[1] - 1, fs = f[0] === "-" ? revComp(seq) : seq;
      for (var a = 0; a < cs.length; a++) {
        if (cs[a].c !== "ATG") continue;
        var z = a + 1; while (z < cs.length && cs[z].a !== "*") z++;
        var stop = z < cs.length, last = stop ? z : cs.length - 1, aa = (stop ? z : cs.length) - a;
        if (aa < minAA) continue;
        var s0 = o + 3 * a, e0 = o + 3 * (last + 1), fw = f[0] === "+";
        res.push({ frame: f, strand: f[0], from: fw ? s0 + 1 : L - e0 + 1, to: fw ? e0 : L - s0, aa: aa, nt: e0 - s0,
          startCodon: true, stopCodon: stop, complete: stop, dna: fs.slice(s0, e0),
          protein: cs.slice(a, last + 1).map(function (x) { return x.a; }).join(""), ci0: a, ci1: last });
      }
    });
    return res.sort(function (x, y) { return y.aa - x.aa || x.from - y.from; });
  }

  // GC% = 100*(G+C)/(A+C+G+T). Ambiguous bases are rejected upstream, so the denominator equals the length.
  function stats(seq) {
    var c = { A: 0, C: 0, G: 0, T: 0 };
    for (var n = 0; n < seq.length; n++) if (c[seq[n]] !== undefined) c[seq[n]]++;
    var d = c.A + c.C + c.G + c.T;
    return { length: seq.length, counts: c, gc: d ? 100 * (c.G + c.C) / d : 0, codons: Math.floor(seq.length / 3) };
  }
  function aaComposition(protein) { // stop codons excluded; length = residues
    var m = {}, n = 0;
    for (var q = 0; q < protein.length; q++) if (protein[q] !== "*") { m[protein[q]] = (m[protein[q]] || 0) + 1; n++; }
    return { counts: m, length: n };
  }
  return { FRAMES: FRAMES, CODON_TABLE: CT, parse: parse, revComp: revComp, transcribe: transcribe, translateFrame: translateFrame, findOrfs: findOrfs, stats: stats, aaComposition: aaComposition };
});
