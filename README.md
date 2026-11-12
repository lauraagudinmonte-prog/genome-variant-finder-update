# genome-variant-finder

[![tests](https://github.com/YOUR-USERNAME/genome-variant-finder/actions/workflows/tests.yml/badge.svg)](https://github.com/YOUR-USERNAME/genome-variant-finder/actions/workflows/tests.yml)

A small command-line toolkit for DNA sequences with two tools:

1. **Genome comparison** (`genome_compare`): compares two FASTA sequences position by position and lists the single-base differences, classified as transitions or transversions.
2. **Motif search** (`genome_compare.motif`): finds a motif such as `TATANNNATGNNNT` on both the forward and the reverse-complement strand.

It started from coursework exercises in a Bioinformatics master's programme (comparing SARS-CoV-2 genomes and searching for sequence patterns). I rewrote them as tested, reusable tools with no external dependencies (Python 3.10+ standard library only).

## Tool 1: genome comparison

### Usage

```bash
python -m genome_compare REFERENCE.fasta SAMPLE.fasta [--ignore-ambiguous] [-o differences.tsv]
```

| Option | Meaning |
| --- | --- |
| `--ignore-ambiguous` | Skip positions where either base is not A, C, G or T (for example `N`) instead of reporting them. |
| `-o, --output` | Also write the differences to a tab-separated file. |

#### Try it with the included demo data

The files in `examples/` are **synthetic** sequences (not real genomes) so the tool can be run without downloading anything:

```bash
python -m genome_compare examples/reference_demo.fasta examples/sample_demo.fasta
```

```text
Reference: DEMO_REF (150 bp)
Sample:    DEMO_SAMPLE (140 bp)
Compared:  first 140 positions
Note:      lengths differ; positions beyond the shorter one were not compared
Differences: 4 (2.8571% of 140 evaluated positions)
  transitions: 2, transversions: 1, other: 1

position	reference	alternate	type
20	G	A	transition
75	A	G	transition
110	T	G	transversion
130	A	N	other
```

With `--ignore-ambiguous`, the `N` at position 130 is skipped and counted separately (3 differences over 139 evaluated positions).

#### Using real genomes

Download two sequences in FASTA format from [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore/), for example the SARS-CoV-2 reference `NC_045512.2` and a sample such as `MN985325.1`, save them in a local `data/` folder (ignored by git), and run:

```bash
python -m genome_compare data/NC_045512.2.fasta data/MN985325.1.fasta -o data/differences.tsv
```

#### How it works

1. `read_fasta` reads the first record of each file, removes line breaks and whitespace, and upper-cases the sequence.
2. `compare_sequences` walks both sequences with `zip` (so it stops at the shorter one) and records every position where the bases differ.
3. Each difference is classified: **transition** (A↔G or C↔T), **transversion** (purine↔pyrimidine) or **other** (ambiguity codes such as `N`).

Positions are **1-based**, as in genome browsers and NCBI.

#### Limitations

- This is a position-by-position comparison, **not an alignment**. An insertion or deletion shifts every following position, so real indels are reported as many spurious differences. For divergent sequences, align first (e.g. MAFFT or Biopython's pairwise aligner).
- Only the first record of each FASTA file is used.
- Only nucleotide sequences are supported.

## Tool 2: motif search on both strands

Searches a sequence for a motif where `N` means "any one of A, C, G, T". The default motif is `TATANNNATGNNNT` (TATA, any three bases, ATG, any three bases, T). The sequence is searched on the forward strand and on its reverse complement, and every hit is reported in **forward-strand coordinates** (1-based, inclusive).

```bash
python -m genome_compare.motif GGCCTATACCCATGGGGTAACC
python -m genome_compare.motif --fasta data/genome.fasta --motif GANNTC
```

```text
Motif TATANNNATGNNNT found 1 time(s):
  forward strand: positions 5-18 (forward coordinates), reads TATACCCATGGGGT
```

`python examples/motif_examples.py` runs three demonstration cases: a motif on the forward strand, one on the reverse strand, and a sequence with no match on either.

How it works:

- A regular expression is built from the motif, with `N` turned into `[ACGT]`. It uses a lookahead, so **overlapping** occurrences are found as well.
- The same search is repeated on the reverse complement. A hit at index `i` of the reverse complement maps back to forward positions `len - i - motif_length + 1` to `len - i`.
- `N` in the *sequence* is never treated as a wildcard match: the motif needs real A, C, G or T bases there. The reverse complement keeps `N` as `N`, because real genomes contain them.
- Invalid characters in the sequence or the motif raise a clear error.

## Tests

```bash
python -m unittest discover -s tests -v
```

The tests cover both tools: 1-based positions, case handling, unequal lengths, ambiguous bases, transition/transversion classification, FASTA parsing errors, reverse-strand coordinate mapping, overlapping motif hits, input validation and the command-line interfaces.

## Possible next steps

- Pairwise alignment before comparing, to handle indels.
- Multi-record FASTA input and comparison of several samples against one reference.
- Report variants in VCF format.
- Support IUPAC ambiguity codes in motifs (R, Y, W...).

## License

MIT, see [LICENSE](LICENSE).
