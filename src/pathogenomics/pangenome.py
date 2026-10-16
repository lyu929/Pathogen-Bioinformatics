"""Pangenome analysis from Roary's ``gene_presence_absence.csv``."""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

#: Roary's metadata columns; every column after these is one isolate
ROARY_METADATA = [
    "Gene",
    "Non-unique Gene name",
    "Annotation",
    "No. isolates",
    "No. sequences",
    "Avg sequences per isolate",
    "Genome Fragment",
    "Order within Fragment",
    "Accessory Fragment",
    "Accessory Order with Fragment",
    "QC",
    "Min group size nuc",
    "Max group size nuc",
    "Avg group size nuc",
]

#: Roary's default partition boundaries (fraction of isolates carrying a gene cluster)
PARTITIONS = (("core", 0.99), ("soft_core", 0.95), ("shell", 0.15), ("cloud", 0.0))


def read_roary(path: str | Path) -> pd.DataFrame:
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = [c for c in ROARY_METADATA if c not in table.columns]
    if missing:
        raise ValueError(f"{path}: not a Roary gene_presence_absence.csv (missing {missing})")
    table["No. isolates"] = table["No. isolates"].astype(int)
    return table.set_index("Gene")


def isolate_columns(table: pd.DataFrame) -> list[str]:
    meta = set(ROARY_METADATA)
    return [c for c in table.columns if c not in meta]


def infer_isolate_names(table: pd.DataFrame) -> dict[str, str]:
    """Map Galaxy column names (``Prokka_on_dataset_102__gff``) to the locus-tag prefix
    Prokka used for that genome (``EN93``), which is how the isolates are known."""
    names: dict[str, str] = {}
    for col in isolate_columns(table):
        prefixes = Counter(
            re.split(r"_\d", v.split("\t")[0])[0]
            for v in table[col]
            if v and re.match(r"^[A-Za-z0-9]+_\d", v)
        )
        names[col] = prefixes.most_common(1)[0][0] if prefixes else col
    if len(set(names.values())) != len(names):
        raise ValueError(f"locus-tag prefixes are not unique: {names}")
    return names


def presence_matrix(table: pd.DataFrame, rename: bool = True) -> pd.DataFrame:
    """Boolean gene-cluster x isolate matrix (a non-empty cell means present)."""
    cols = isolate_columns(table)
    presence = table[cols].ne("")
    if rename:
        presence = presence.rename(columns=infer_isolate_names(table))
    if not (presence.sum(axis=1) == table["No. isolates"]).all():
        raise ValueError("presence matrix disagrees with Roary's 'No. isolates' column")
    return presence


def partition(presence: pd.DataFrame, boundaries=PARTITIONS) -> pd.Series:
    """Assign each gene cluster to core / soft_core / shell / cloud by isolate fraction."""
    frac = presence.sum(axis=1) / presence.shape[1]
    labels = pd.Series(index=presence.index, dtype=object)
    remaining = pd.Series(True, index=presence.index)
    for name, lower in boundaries:
        hit = remaining & (frac >= lower - 1e-12)
        labels[hit] = name
        remaining &= ~hit
    return labels


def partition_summary(presence: pd.DataFrame) -> pd.DataFrame:
    labels = partition(presence)
    order = [name for name, _ in PARTITIONS]
    counts = labels.value_counts().reindex(order, fill_value=0)
    counts.loc["total"] = counts.sum()
    out = pd.DataFrame(
        {"clusters": counts.astype(int), "pct": (100 * counts / counts["total"]).round(2)}
    )
    out.index.name = "partition"
    return out


def rarefaction(presence: pd.DataFrame, permutations: int = 200, seed: int = 0) -> pd.DataFrame:
    """Pan- and core-genome size as genomes are added in random order.

    Returns mean and standard deviation over ``permutations`` random orderings for
    k = 1..N genomes.
    """
    rng = np.random.default_rng(seed)
    mat = presence.to_numpy(dtype=bool)
    n = mat.shape[1]
    pan = np.zeros((permutations, n), dtype=np.int64)
    core = np.zeros((permutations, n), dtype=np.int64)
    for p in range(permutations):
        order = rng.permutation(n)
        seen_any = np.zeros(mat.shape[0], dtype=bool)
        seen_all = np.ones(mat.shape[0], dtype=bool)
        for k, col in enumerate(order):
            seen_any |= mat[:, col]
            seen_all &= mat[:, col]
            pan[p, k] = seen_any.sum()
            core[p, k] = seen_all.sum()
    return pd.DataFrame(
        {
            "genomes": np.arange(1, n + 1),
            "pan_mean": pan.mean(axis=0),
            "pan_sd": pan.std(axis=0),
            "core_mean": core.mean(axis=0),
            "core_sd": core.std(axis=0),
        }
    )


def jaccard_distance(presence: pd.DataFrame, accessory_only: bool = True) -> pd.DataFrame:
    """Pairwise Jaccard distance between isolates on (accessory) gene content."""
    mat = presence
    if accessory_only:
        mat = presence[~presence.all(axis=1)]
    names = list(mat.columns)
    arr = mat.to_numpy(dtype=bool)
    out = np.zeros((len(names), len(names)))
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            union = np.logical_or(arr[:, i], arr[:, j]).sum()
            inter = np.logical_and(arr[:, i], arr[:, j]).sum()
            out[i, j] = out[j, i] = 1 - inter / union if union else 0.0
    return pd.DataFrame(out, index=names, columns=names)


def unique_genes(table: pd.DataFrame, presence: pd.DataFrame) -> pd.DataFrame:
    """Gene clusters found in exactly one isolate, with their Roary annotation."""
    single = presence[presence.sum(axis=1) == 1]
    owner = single.idxmax(axis=1)
    return pd.DataFrame(
        {"isolate": owner, "annotation": table.loc[single.index, "Annotation"]}
    ).sort_values(["isolate", "annotation"])


def absent_in(presence: pd.DataFrame, table: pd.DataFrame, isolate: str) -> pd.DataFrame:
    """Clusters present in every isolate except ``isolate`` (e.g. genes lost relative to a reference)."""
    others = [c for c in presence.columns if c != isolate]
    mask = presence[others].all(axis=1) & ~presence[isolate]
    return table.loc[mask[mask].index, ["Annotation"]]
