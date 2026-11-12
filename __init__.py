"""genome_compare: small tools for DNA sequences.

Sequence comparison lives in ``genome_compare.core`` (CLI: ``python -m genome_compare``).
Motif search lives in ``genome_compare.motif`` (CLI: ``python -m genome_compare.motif``).
"""

from .core import Comparison, FastaRecord, Variant, compare_sequences, read_fasta

__all__ = ["Comparison", "FastaRecord", "Variant", "compare_sequences", "read_fasta"]
__version__ = "0.2.0"
