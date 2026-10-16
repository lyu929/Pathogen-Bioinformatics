"""Bulk RNA-seq: count matrix QC, normalisation, PCA and differential expression.

Differential expression is run with PyDESeq2 (a Python port of DESeq2) and then
*cross-checked* against the DESeq2 results produced in Galaxy for the course, so
the two independent implementations validate each other.
"""

from __future__ import annotations

import contextlib
import io
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from pathogenomics.io.galaxy import DESEQ2_COLUMNS, read_featurecounts_dir


@dataclass
class RnaSeqDataset:
    name: str
    counts: pd.DataFrame  # genes x samples (int)
    samples: pd.DataFrame  # indexed by sample, must contain the design factor

    def __post_init__(self) -> None:
        missing = set(self.samples.index) ^ set(self.counts.columns)
        if missing:
            raise ValueError(
                f"{self.name}: samples sheet and count matrix disagree on {sorted(missing)}"
            )
        self.counts = self.counts[list(self.samples.index)]


def load_dataset(name: str, counts_dir: str | Path, samples_tsv: str | Path) -> RnaSeqDataset:
    counts = read_featurecounts_dir(counts_dir, "*.tsv")
    samples = pd.read_csv(samples_tsv, sep="\t", dtype=str).set_index("sample")
    return RnaSeqDataset(name=name, counts=counts, samples=samples)


def median_of_ratios_size_factors(counts: pd.DataFrame) -> pd.Series:
    """DESeq2 size factors: median ratio of each sample to the per-gene geometric mean.

    Only genes with a non-zero count in every sample enter the geometric mean
    (``estimateSizeFactors`` default).
    """
    arr = counts.to_numpy(dtype=float)
    expressed = (arr > 0).all(axis=1)
    if not expressed.any():
        raise ValueError("no gene is expressed in every sample; median-of-ratios is undefined")
    log_counts = np.log(arr[expressed])
    log_geo_mean = log_counts.mean(axis=1, keepdims=True)
    factors = np.exp(np.median(log_counts - log_geo_mean, axis=0))
    return pd.Series(factors, index=counts.columns, name="size_factor")


def normalized_counts(counts: pd.DataFrame, size_factors: pd.Series | None = None) -> pd.DataFrame:
    sf = median_of_ratios_size_factors(counts) if size_factors is None else size_factors
    return counts.div(sf, axis=1)


def library_summary(dataset: RnaSeqDataset, min_count: int = 10) -> pd.DataFrame:
    counts = dataset.counts
    total = counts.sum()
    top = counts.apply(lambda c: c.nlargest(1).iloc[0]) / total
    out = pd.DataFrame(
        {
            "assigned_reads": total,
            "genes_detected": (counts >= min_count).sum(),
            "top_gene_pct": (100 * top).round(2),
            "size_factor": median_of_ratios_size_factors(counts).round(4),
        }
    )
    return dataset.samples.join(out)


def log_normalized(counts: pd.DataFrame, pseudocount: float = 1.0) -> pd.DataFrame:
    return np.log2(normalized_counts(counts) + pseudocount)


def pca(
    log_expr: pd.DataFrame, n_top: int = 500, n_components: int = 2
) -> tuple[pd.DataFrame, np.ndarray]:
    """PCA of samples on the ``n_top`` most variable genes (as ``DESeq2::plotPCA``).

    Returns ``(coordinates, explained_variance_ratio)``.
    """
    variances = log_expr.var(axis=1)
    top = log_expr.loc[variances.nlargest(min(n_top, len(variances))).index]
    centered = top.sub(top.mean(axis=1), axis=0).to_numpy().T  # samples x genes
    u, s, _ = np.linalg.svd(centered, full_matrices=False)
    coords = u[:, :n_components] * s[:n_components]
    # fix the arbitrary SVD sign so results are deterministic across BLAS builds
    signs = np.sign(coords[np.abs(coords).argmax(axis=0), np.arange(coords.shape[1])])
    coords = coords * signs
    explained = (s**2) / (s**2).sum()
    cols = [f"PC{i + 1}" for i in range(n_components)]
    return pd.DataFrame(coords, index=log_expr.columns, columns=cols), explained[:n_components]


def run_deseq2(
    dataset: RnaSeqDataset,
    factor: str,
    reference: str,
    test: str,
    alpha: float = 0.05,
    n_cpus: int = 1,
) -> tuple[pd.DataFrame, pd.Series]:
    """Wald test of ``test`` vs ``reference`` for ``~factor`` with PyDESeq2.

    Returns ``(results, size_factors)``; ``results`` has DESeq2's column names.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    logging.getLogger("pydeseq2").setLevel(logging.ERROR)
    inference = DefaultInference(n_cpus=n_cpus)
    metadata = dataset.samples[[factor]].copy()
    levels = set(metadata[factor])
    if {reference, test} - levels:
        raise ValueError(
            f"levels {reference!r}/{test!r} not in design factor {factor!r} ({sorted(levels)})"
        )
    dds = DeseqDataSet(
        counts=dataset.counts.T,
        metadata=metadata,
        design=f"~{factor}",
        refit_cooks=True,
        inference=inference,
        quiet=True,
    )
    with contextlib.redirect_stdout(io.StringIO()):
        dds.deseq2()
        result = DeseqStats(
            dds, contrast=[factor, test, reference], alpha=alpha, inference=inference, quiet=True
        )
        result.summary()
    table = result.results_df[DESEQ2_COLUMNS].copy()
    table.index.name = "gene_id"
    size_factors = pd.Series(
        np.asarray(dds.obs["size_factors"]), index=dataset.counts.columns, name="size_factor"
    )
    return table, size_factors


def significant(
    results: pd.DataFrame, alpha: float = 0.05, min_abs_lfc: float = 1.0
) -> pd.DataFrame:
    mask = (results["padj"] < alpha) & (results["log2FoldChange"].abs() >= min_abs_lfc)
    return results[mask.fillna(False)]


def compare_results(
    ours: pd.DataFrame, reference: pd.DataFrame, alpha: float = 0.05, min_abs_lfc: float = 1.0
) -> dict[str, float | int]:
    """Concordance of two DESeq2 result tables (e.g. PyDESeq2 vs Galaxy/R DESeq2)."""
    common = ours.index.intersection(reference.index)
    a, b = ours.loc[common], reference.loc[common]
    finite = a["log2FoldChange"].notna() & b["log2FoldChange"].notna()
    pa = -np.log10(a["pvalue"].clip(lower=1e-300))
    pb = -np.log10(b["pvalue"].clip(lower=1e-300))
    pfin = pa.notna() & pb.notna()
    sig_a = set(significant(a, alpha, min_abs_lfc).index)
    sig_b = set(significant(b, alpha, min_abs_lfc).index)
    both = sig_a & sig_b
    union = sig_a | sig_b
    same_sign = sum(
        1 for g in both if np.sign(a.at[g, "log2FoldChange"]) == np.sign(b.at[g, "log2FoldChange"])
    )
    tested = a["padj"].notna() & b["padj"].notna()
    big = tested & (a["log2FoldChange"].abs() > min_abs_lfc)
    return {
        "genes_compared": int(len(common)),
        "basemean_pearson": round(float(stats.pearsonr(a["baseMean"], b["baseMean"])[0]), 6),
        "lfc_pearson": round(
            float(
                stats.pearsonr(a.loc[finite, "log2FoldChange"], b.loc[finite, "log2FoldChange"])[0]
            ),
            5,
        ),
        "lfc_max_abs_diff": round(
            float((a.loc[finite, "log2FoldChange"] - b.loc[finite, "log2FoldChange"]).abs().max()),
            4,
        ),
        "lfc_pearson_tested_genes": round(
            float(
                stats.pearsonr(a.loc[tested, "log2FoldChange"], b.loc[tested, "log2FoldChange"])[0]
            ),
            5,
        ),
        "lfc_magnitude_ratio_ref_over_ours": round(
            float(
                (b.loc[big, "log2FoldChange"].abs() / a.loc[big, "log2FoldChange"].abs()).median()
            ),
            3,
        )
        if big.any()
        else float("nan"),
        "neglog10p_spearman": round(float(stats.spearmanr(pa[pfin], pb[pfin])[0]), 5),
        "padj_na_ours": int(a["padj"].isna().sum()),
        "padj_na_reference": int(b["padj"].isna().sum()),
        "significant_ours": len(sig_a),
        "significant_reference": len(sig_b),
        "significant_both": len(both),
        "jaccard": round(len(both) / len(union), 4) if union else 1.0,
        "sign_agreement_pct": round(100 * same_sign / len(both), 2) if both else 100.0,
    }


def annotate(results: pd.DataFrame, annotation: pd.DataFrame) -> pd.DataFrame:
    """Join gene symbol and product from :func:`pathogenomics.io.gff.gene_annotation_table`."""
    cols = [c for c in ("gene", "product") if c in annotation.columns]
    return results.join(annotation[cols], how="left")


def direction_counts(
    results: pd.DataFrame, alpha: float = 0.05, min_abs_lfc: float = 1.0
) -> dict[str, int]:
    sig = significant(results, alpha, min_abs_lfc)
    return {
        "up": int((sig["log2FoldChange"] > 0).sum()),
        "down": int((sig["log2FoldChange"] < 0).sum()),
    }
