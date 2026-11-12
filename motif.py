"""Search a DNA sequence for a motif on both strands.

The default motif is ``TATANNNATGNNNT``: TATA, any three bases, ATG, any three
bases, T. In a motif, ``N`` stands for any one of A, C, G, T.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from .core import read_fasta

DEFAULT_MOTIF = "TATANNNATGNNNT"

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


@dataclass(frozen=True)
class MotifHit:
    """One occurrence of a motif.

    ``start`` and ``end`` are 1-based, inclusive coordinates on the *forward*
    strand of the sequence that was searched, for both strands. ``match`` is
    the motif as it reads 5'->3' on the strand where it was found.
    """

    strand: str  # '+' (forward) or '-' (reverse complement)
    start: int
    end: int
    match: str


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence.

    ``N`` (unknown base) is kept as ``N``, since real genomes contain them.

    Raises:
        ValueError: if the sequence contains anything other than A, C, G, T, N.
    """
    sequence = sequence.upper()
    invalid = set(sequence) - set("ACGTN")
    if invalid:
        raise ValueError(f"invalid DNA characters: {', '.join(sorted(invalid))}")
    return sequence.translate(_COMPLEMENT)[::-1]


def _motif_to_regex(motif: str) -> re.Pattern[str]:
    motif = motif.upper()
    invalid = set(motif) - set("ACGTN")
    if not motif or invalid:
        raise ValueError("motif must be non-empty and contain only A, C, G, T and N")
    body = "".join("[ACGT]" if base == "N" else base for base in motif)
    # The lookahead lets overlapping occurrences be found too.
    return re.compile(f"(?=({body}))")


def find_motif(sequence: str, motif: str = DEFAULT_MOTIF) -> list[MotifHit]:
    """Find every occurrence of ``motif`` on the forward and reverse strands.

    Hits are returned sorted by start position. The motif is searched on the
    reverse complement of ``sequence`` as well, and those hits are converted
    back to forward-strand coordinates.

    Raises:
        ValueError: if the sequence or the motif contains invalid characters.
    """
    sequence = sequence.upper()
    pattern = _motif_to_regex(motif)
    length = len(motif)
    total = len(sequence)

    hits = [
        MotifHit("+", m.start() + 1, m.start() + length, m.group(1))
        for m in pattern.finditer(sequence)
    ]
    # Search the reverse strand; index i on it maps to forward coordinates
    # total - i - length + 1 ... total - i.
    for m in pattern.finditer(reverse_complement(sequence)):
        i = m.start()
        hits.append(MotifHit("-", total - i - length + 1, total - i, m.group(1)))

    return sorted(hits, key=lambda hit: (hit.start, hit.strand))


def describe(sequence: str, motif: str = DEFAULT_MOTIF) -> str:
    """Human-readable report of where the motif occurs in a sequence."""
    hits = find_motif(sequence, motif)
    if not hits:
        return f"Motif {motif} not found on either strand."
    names = {"+": "forward strand", "-": "reverse strand"}
    lines = [f"Motif {motif} found {len(hits)} time(s):"]
    for hit in hits:
        lines.append(
            f"  {names[hit.strand]}: positions {hit.start}-{hit.end} "
            f"(forward coordinates), reads {hit.match}"
        )
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="genome_compare.motif",
        description="Search a DNA sequence for a motif on both strands.",
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("sequence", nargs="?", help="DNA sequence given directly")
    source.add_argument("--fasta", type=Path, help="FASTA file (first record is used)")
    parser.add_argument(
        "--motif",
        default=DEFAULT_MOTIF,
        help=f"motif using A, C, G, T and N for any base (default: {DEFAULT_MOTIF})",
    )
    args = parser.parse_args(argv)

    try:
        sequence = read_fasta(args.fasta).sequence if args.fasta else args.sequence
        print(describe(sequence, args.motif))
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
