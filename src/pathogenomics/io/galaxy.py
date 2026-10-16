"""Readers for Galaxy tabular outputs (featureCounts, DESeq2)."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

DESEQ2_COLUMNS = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]


def read_featurecounts(path: str | Path, sample: str | None = None) -> pd.Series:
    """Read a Galaxy featureCounts ``Counts`` table (``Geneid<TAB>count``) as an int Series.

    Galaxy names the count column after the BAM dataset (e.g. ``"71 Bowtie2 on dataset
    85, 86, and 103: alignments"``); the series is renamed to ``sample`` (default: the
    first ``SRX/SRR/ERR`` accession in the file name, else the file stem).
    """
    path = Path(path)
    table = pd.read_csv(path, sep="\t", comment="#", dtype={0: str})
    if table.shape[1] == 2:
        gene_col, count_col = table.columns
    elif {"Geneid", "Chr", "Start", "End", "Strand", "Length"}.issubset(table.columns):
        # full featureCounts layout: counts are the last column
        gene_col, count_col = "Geneid", table.columns[-1]
    else:
        raise ValueError(f"{path}: not a featureCounts table (columns={list(table.columns)})")
    if sample is None:
        match = re.search(r"((?:SRX|SRR|ERR|DRR)\d+)", path.name)
        sample = match.group(1) if match else path.name.split(".")[0]
    counts = table.set_index(gene_col)[count_col]
    if counts.index.duplicated().any():
        raise ValueError(f"{path}: duplicated gene identifiers")
    if (counts < 0).any():
        raise ValueError(f"{path}: negative counts")
    counts = counts.astype("int64")
    counts.index.name = "gene_id"
    counts.name = sample
    return counts


def read_featurecounts_dir(directory: str | Path, pattern: str = "*.tsv") -> pd.DataFrame:
    """Combine per-sample featureCounts tables into a genes x samples matrix."""
    files = sorted(Path(directory).glob(pattern))
    if not files:
        raise FileNotFoundError(f"no featureCounts tables matching {pattern!r} in {directory}")
    series = [read_featurecounts(f) for f in files]
    index = series[0].index
    for s in series[1:]:
        if not s.index.equals(index):
            raise ValueError(f"gene set/order of {s.name} differs from {series[0].name}")
    return pd.concat(series, axis=1)


def read_deseq2_table(path: str | Path) -> pd.DataFrame:
    """Read a DESeq2 results table (with or without Galaxy's missing header)."""
    path = Path(path)
    first = path.open().readline().split("\t")
    has_header = any(c.strip() in DESEQ2_COLUMNS for c in first)
    table = pd.read_csv(path, sep="\t", header=0 if has_header else None)
    if not has_header:
        table.columns = ["gene_id", *DESEQ2_COLUMNS][: table.shape[1]]
    table = table.rename(columns={table.columns[0]: "gene_id"}).set_index("gene_id")
    missing = set(DESEQ2_COLUMNS) - set(table.columns)
    if missing:
        raise ValueError(f"{path}: missing DESeq2 columns {sorted(missing)}")
    return table[DESEQ2_COLUMNS]
