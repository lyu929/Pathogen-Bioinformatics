"""De novo assembly statistics (QUAST-style N50/L50/NG50) for SPAdes/Shovill contigs."""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from pathogenomics.io.seq import read_fasta

_SPADES_COV = re.compile(r"_cov_([0-9.]+)")
_SHOVILL_COV = re.compile(r"\bcov=([0-9.]+)")


def nx(
    lengths: Sequence[int], fraction: float = 0.5, genome_size: int | None = None
) -> tuple[int, int]:
    """Return ``(Nx, Lx)`` for ``fraction`` in (0, 1].

    ``Nx`` is the length of the shortest contig in the smallest set of longest contigs
    covering ``fraction`` of the assembly (or of ``genome_size`` for NGx); ``Lx`` is the
    number of contigs in that set. Returns ``(0, 0)`` if the target is never reached.
    """
    if not 0 < fraction <= 1:
        raise ValueError("fraction must be in (0, 1]")
    ordered = np.sort(np.asarray(lengths, dtype=np.int64))[::-1]
    if ordered.size == 0:
        return 0, 0
    target = fraction * (genome_size if genome_size is not None else int(ordered.sum()))
    cumulative = np.cumsum(ordered)
    idx = int(np.searchsorted(cumulative, target, side="left"))
    if idx >= ordered.size:
        return 0, 0
    return int(ordered[idx]), idx + 1


def contig_coverage(header: str) -> float | None:
    """k-mer/read coverage encoded in SPAdes (``_cov_24.8``) or Shovill (``cov=58.2``) headers."""
    for pattern in (_SPADES_COV, _SHOVILL_COV):
        match = pattern.search(header)
        if match:
            return float(match.group(1))
    return None


@dataclass(frozen=True)
class AssemblyStats:
    name: str
    contigs: int
    contigs_ge_1kb: int
    total_length: int
    largest_contig: int
    n50: int
    l50: int
    n90: int
    l90: int
    gc_percent: float
    n_per_100kb: float
    mean_coverage: float | None
    ng50: int | None = None

    def as_dict(self) -> dict[str, object]:
        return self.__dict__.copy()


def assembly_stats(
    path: str | Path, name: str | None = None, min_length: int = 0, genome_size: int | None = None
) -> AssemblyStats:
    """Compute assembly metrics from a FASTA file, ignoring contigs shorter than ``min_length``.

    ``mean_coverage`` is the length-weighted mean of the coverage values found in the
    contig headers (``None`` if headers carry no coverage).
    """
    lengths: list[int] = []
    gc = acgt = n_count = 0
    cov_weighted = 0.0
    cov_length = 0
    for header, seq in read_fasta(path):
        if len(seq) < min_length:
            continue
        s = seq.upper()
        lengths.append(len(s))
        g_c = s.count("G") + s.count("C")
        gc += g_c
        acgt += g_c + s.count("A") + s.count("T")
        n_count += s.count("N")
        cov = contig_coverage(header)
        if cov is not None:
            cov_weighted += cov * len(s)
            cov_length += len(s)
    if not lengths:
        raise ValueError(f"{path}: no contigs >= {min_length} bp")
    total = int(sum(lengths))
    n50, l50 = nx(lengths, 0.5)
    n90, l90 = nx(lengths, 0.9)
    ng50 = nx(lengths, 0.5, genome_size)[0] if genome_size else None
    return AssemblyStats(
        name=name or Path(path).name.split(".")[0],
        contigs=len(lengths),
        contigs_ge_1kb=sum(1 for x in lengths if x >= 1000),
        total_length=total,
        largest_contig=max(lengths),
        n50=n50,
        l50=l50,
        n90=n90,
        l90=l90,
        gc_percent=round(100.0 * gc / acgt, 2) if acgt else 0.0,
        n_per_100kb=round(100_000 * n_count / total, 2),
        mean_coverage=round(cov_weighted / cov_length, 1) if cov_length else None,
        ng50=ng50,
    )


def assembly_table(
    paths: Sequence[str | Path], min_length: int = 0, genome_size: int | None = None
) -> pd.DataFrame:
    rows = [
        assembly_stats(p, min_length=min_length, genome_size=genome_size).as_dict() for p in paths
    ]
    return pd.DataFrame(rows).set_index("name")


def contig_lengths(path: str | Path, min_length: int = 0) -> np.ndarray:
    return np.asarray([len(s) for _, s in read_fasta(path) if len(s) >= min_length], dtype=np.int64)
