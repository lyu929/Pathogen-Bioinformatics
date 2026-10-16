"""Read quality control: FASTQ statistics and Cutadapt-compatible quality trimming.

``quality_trim_index`` implements the modified BWA algorithm used by Cutadapt's
``--quality-cutoff`` option, so that trimming the course FASTQ with this module
reproduces the Galaxy Cutadapt report base-for-base (see ``tests/test_qc.py``).
"""

from __future__ import annotations

import contextlib
import gzip
from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from pathogenomics.io.reports import CutadaptReport, FastQCReport
from pathogenomics.io.seq import read_fastq


def quality_trim_index(
    qualities: str, cutoff_back: int, cutoff_front: int = 0, base: int = 33
) -> tuple[int, int]:
    """Return ``(start, stop)`` so that ``read[start:stop]`` is the quality-trimmed read.

    For the 3' end, the trim point is the position that maximises
    ``sum(cutoff - q_i)`` over the read suffix, scanning from the end and stopping
    as soon as the running sum becomes negative (Cutadapt/BWA). The 5' end is
    handled symmetrically with ``cutoff_front`` (0 disables it).
    """
    n = len(qualities)
    start, stop = 0, n
    if cutoff_front > 0:
        running = best = 0
        for i in range(n):
            running += cutoff_front - (ord(qualities[i]) - base)
            if running < 0:
                break
            if running > best:
                best, start = running, i + 1
    if cutoff_back > 0:
        running = best = 0
        for i in range(n - 1, -1, -1):
            running += cutoff_back - (ord(qualities[i]) - base)
            if running < 0:
                break
            if running > best:
                best, stop = running, i
    if start >= stop:
        return 0, 0
    return start, stop


@dataclass(frozen=True)
class TrimSummary:
    reads_processed: int
    too_short: int
    written: int
    bp_processed: int
    bp_quality_trimmed: int
    bp_written: int

    def matches(self, report: CutadaptReport) -> bool:
        """True if every count equals the corresponding Cutadapt report field."""
        return all(getattr(self, k) == getattr(report, k) for k in asdict(self))


def trim_fastq(
    path: str | Path,
    quality_cutoff: int = 20,
    min_length: int = 22,
    front_cutoff: int = 0,
    output: str | Path | None = None,
    base: int = 33,
) -> TrimSummary:
    """Quality-trim a single-end FASTQ and drop reads shorter than ``min_length``.

    Mirrors ``cutadapt --quality-cutoff=[front,]back --minimum-length=N`` without
    adapter removal. If ``output`` is given, passing reads are written there
    (gzip-compressed when the name ends in ``.gz``).
    """
    reads = short = written = bp = bp_trim = bp_out = 0
    with contextlib.ExitStack() as stack:
        out = None
        if output is not None:
            output = Path(output)
            opener = gzip.open if output.suffix == ".gz" else open
            out = stack.enter_context(opener(output, "wt"))
        for rec in read_fastq(path):
            reads += 1
            bp += len(rec.sequence)
            start, stop = quality_trim_index(rec.quality, quality_cutoff, front_cutoff, base)
            kept = stop - start
            bp_trim += len(rec.sequence) - kept
            if kept < min_length:
                short += 1
                continue
            written += 1
            bp_out += kept
            if out is not None:
                out.write(
                    f"@{rec.name}\n{rec.sequence[start:stop]}\n+\n{rec.quality[start:stop]}\n"
                )
    return TrimSummary(reads, short, written, bp, bp_trim, bp_out)


@dataclass
class FastqStats:
    reads: int
    bases: int
    min_length: int
    max_length: int
    gc_percent: float
    n_bases: int
    per_position_mean_quality: np.ndarray
    per_position_count: np.ndarray
    mean_read_quality: np.ndarray

    def summary(self) -> dict[str, float | int]:
        return {
            "reads": self.reads,
            "bases": self.bases,
            "min_length": self.min_length,
            "max_length": self.max_length,
            "gc_percent": round(self.gc_percent, 2),
            "n_bases": self.n_bases,
            "mean_quality": round(float(self.mean_read_quality.mean()), 2),
            "pct_reads_meanQ_ge_20": round(float((self.mean_read_quality >= 20).mean() * 100), 2),
        }


def fastq_stats(path: str | Path, base: int = 33) -> FastqStats:
    """Single-pass FASTQ summary comparable to FastQC's basic and per-base modules."""
    qual_sum = np.zeros(0, dtype=np.int64)
    qual_n = np.zeros(0, dtype=np.int64)
    read_means: list[float] = []
    gc = at = n_bases = reads = bases = 0
    lmin, lmax = 10**9, 0
    for rec in read_fastq(path):
        q = np.frombuffer(rec.quality.encode("ascii"), dtype=np.uint8).astype(np.int64) - base
        length = q.size
        if length > qual_sum.size:
            qual_sum = np.pad(qual_sum, (0, length - qual_sum.size))
            qual_n = np.pad(qual_n, (0, length - qual_n.size))
        qual_sum[:length] += q
        qual_n[:length] += 1
        read_means.append(float(q.mean()) if length else 0.0)
        seq = rec.sequence.upper()
        g_c = seq.count("G") + seq.count("C")
        gc += g_c
        n = seq.count("N")
        n_bases += n
        at += length - g_c - n
        reads += 1
        bases += length
        lmin, lmax = min(lmin, length), max(lmax, length)
    if reads == 0:
        raise ValueError(f"{path}: no reads")
    with np.errstate(invalid="ignore", divide="ignore"):
        per_pos = qual_sum / qual_n
    return FastqStats(
        reads=reads,
        bases=bases,
        min_length=lmin,
        max_length=lmax,
        gc_percent=100.0 * gc / (gc + at) if gc + at else 0.0,
        n_bases=n_bases,
        per_position_mean_quality=per_pos,
        per_position_count=qual_n,
        mean_read_quality=np.asarray(read_means),
    )


def binned_mean_quality(stats: FastqStats, bins: pd.DataFrame) -> np.ndarray:
    """Average quality over FastQC's base groups (``start``/``end`` 1-based), weighting
    each position by the number of reads that reach it (as FastQC pools the bases)."""
    q = stats.per_position_mean_quality
    n = stats.per_position_count
    out = []
    for start, end in zip(bins["start"], bins["end"], strict=True):
        sl = slice(start - 1, end)
        out.append(float(np.average(q[sl], weights=n[sl])))
    return np.asarray(out)


def fastqc_summary(reports: dict[str, FastQCReport]) -> pd.DataFrame:
    """One row per sample: basic statistics plus PASS/WARN/FAIL for each module."""
    rows = []
    for sample, rep in reports.items():
        row: dict[str, object] = {
            "sample": sample,
            "total_sequences": rep.total_sequences,
            "sequence_length": rep.sequence_length,
            "gc_percent": rep.gc_percent,
        }
        row.update(
            {
                module: status.upper()
                for module, status in rep.statuses.items()
                if module != "Basic Statistics"
            }
        )
        rows.append(row)
    return pd.DataFrame(rows).set_index("sample")


def cutadapt_summary(reports: dict[str, CutadaptReport]) -> pd.DataFrame:
    rows = []
    for sample, rep in reports.items():
        rows.append(
            {
                "sample": sample,
                "layout": "paired" if rep.paired else "single",
                "reads_or_pairs": rep.reads_processed,
                "adapter": rep.option("-a") or "",
                "quality_cutoff": rep.option("--quality-cutoff"),
                "min_length": rep.option("--minimum-length"),
                "pct_with_adapter": round(100 * rep.reads_with_adapters / rep.reads_processed, 1)
                if rep.reads_with_adapters is not None
                else None,
                "pct_too_short": round(100 * rep.too_short / rep.reads_processed, 2),
                "pct_reads_kept": round(100 * rep.fraction_reads_kept, 2),
                "pct_bp_quality_trimmed": round(100 * rep.bp_quality_trimmed / rep.bp_processed, 2),
                "pct_bp_kept": round(100 * rep.fraction_bp_kept, 2),
            }
        )
    return pd.DataFrame(rows).set_index("sample")


def flag_failing_modules(summary: pd.DataFrame, modules: Iterable[str] | None = None) -> pd.Series:
    """Comma-separated list of FAIL modules per sample (empty string if none)."""
    cols = [c for c in summary.columns if modules is None or c in set(modules)]
    status_cols = [c for c in cols if summary[c].isin(["PASS", "WARN", "FAIL"]).any()]

    def fails(row: pd.Series) -> str:
        return ", ".join(c for c in status_cols if row.get(c) == "FAIL")

    return summary.apply(fails, axis=1)


def per_base_quality_profile(qualities: Sequence[str], base: int = 33) -> np.ndarray:
    """Mean quality per position for an in-memory list of quality strings (used in tests)."""
    length = max(len(q) for q in qualities)
    total = np.zeros(length)
    count = np.zeros(length)
    for q in qualities:
        arr = np.frombuffer(q.encode("ascii"), dtype=np.uint8) - base
        total[: arr.size] += arr
        count[: arr.size] += 1
    return total / count
