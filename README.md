# DNA Translator

DNA translator and ORF finder for biotechnology. Live tool: https://rehankhanbiotech-pixel.github.io/dna-translator/

## Features
Transcription (mRNA), reverse complement, six-frame translation, ORF finder with minimum-length filter (including nested/overlapping ORFs), GC content and nucleotide counts, amino-acid composition, copy/clear/example controls. Works on mobile and desktop.

## Using the tool
Paste a sequence, press **Translate**. Pick a frame tab to see its protein and codons. Click an ORF row to highlight it and see its DNA and protein.

## Supported input
- DNA (A, C, G, T) in any case; RNA (U is converted to T, with a notice). Mixed T and U is rejected.
- FASTA: header lines (`>`) are ignored; multiline sequences are joined; whitespace and line breaks are removed. If there are several records, only the first is used (with a notice).
- Anything else (digits, `N` and other IUPAC ambiguity codes, symbols) is **rejected** with the character and position. Nothing is silently deleted, so bases on either side of an `N` are never joined.
- Empty input, or a header with no sequence, gives an error message.

## Strand convention
The input is the **coding (sense) strand, 5′→3′**: it has the same sequence as the mRNA (T→U). The **template (antisense) strand** is its reverse complement and is what RNA polymerase reads. Template-strand input is not supported; reverse-complement it first. The "Reverse complement" output is the template strand written 5′→3′.

## Genetic code
Standard genetic code (NCBI table 1) only. Stops are shown as `*`. Start codon for ORFs: `ATG` only (alternative starts and genetic codes are not implemented). Incomplete trailing codons (1–2 nt) are ignored.

## ORF definition and coordinates
- An ORF begins at an `ATG` and ends at the first in-frame stop codon. **Every** ATG is scanned independently, so nested ORFs are all reported (`ATGATGTAA` gives `MM*` and `M*`).
- Complete ORF = start + stop. Incomplete ORF = start but no stop; it runs to the last complete codon and is labelled "no stop".
- Minimum length counts **amino acids before the stop** (stop excluded). Default 30.
- Frames +1/+2/+3 start at base 1/2/3. Frames −1/−2/−3 read the reverse complement from base 1/2/3 of that strand.
- Positions are **1-based, inclusive, on the original input**. The range includes the stop codon. For reverse-strand ORFs, `from < to` in forward coordinates, and the displayed DNA is the reverse-complement (ORF-strand) sequence.

## Statistics
GC% = 100 × (G + C) / (A + C + G + T). Ambiguous bases are rejected, so the denominator equals the sequence length. Empty input never gives NaN/Infinity (GC is 0). Protein length and amino-acid composition exclude stop codons.

## Files
`index.html` (UI), `dna-core.js` (all logic, also used by tests), `test/dna-core.test.js`, `package.json`.

## Running the tests
Requires Node.js 18+ (no dependencies): `npm test` (or `node --test`).

Result of the run on 2026-10-24 (Node v22.22.2): **12 tests, 12 passed, 0 failed.** They cover the cases ATGGCTTAA, ATGATGTAA (nested ORFs), CCCCCC, ATGG, ATGNNNTAA (rejected), AUGGCUUAA (RNA), a reverse-strand ORF with coordinates, minimum-length thresholds, empty/lowercase/FASTA/invalid input, no-stop ORFs, six frames, and 200 random sequences (coordinate round-trip and ORF count check).

Reproducible example: `CCTTAAGCCATG` → one reverse-strand ORF, frame −2, positions 3–11, `ATGGCTTAA`, protein `MA*`.

## Known limitations
No alternative genetic codes or start codons; no ambiguity-code support; circular sequences and ORFs spanning the origin are not considered; only the first FASTA record is processed; the browser UI has not been covered by automated tests (only the logic in `dna-core.js`). This is a teaching/analysis aid, not a validated clinical tool.
