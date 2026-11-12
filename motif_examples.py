"""Three demonstration cases for the TATANNNATGNNNT motif search.

Run from the repository root:  python examples/motif_examples.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from genome_compare.motif import describe, reverse_complement  # noqa: E402

MOTIF = "TATACCCATGGGGT"  # one concrete instance of TATANNNATGNNNT

CASES = {
    "Case 1: motif on the forward strand": "GGCC" + MOTIF + "AACC",
    "Case 2: motif on the reverse-complement strand": "GGCC" + reverse_complement(MOTIF) + "AACC",
    "Case 3: motif on neither strand": "GGCCTATACCCATGCCCAAAAGGTT",
}

for title, sequence in CASES.items():
    print(title)
    print(f"  sequence: {sequence}")
    print("  " + describe(sequence).replace("\n", "\n  "))
    print()
