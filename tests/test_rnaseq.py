import numpy as np
import pandas as pd
import pytest

from pathogenomics.io.galaxy import read_deseq2_table
from pathogenomics.rnaseq import (
    RnaSeqDataset,
    compare_results,
    direction_counts,
    load_dataset,
    log_normalized,
    median_of_ratios_size_factors,
    normalized_counts,
    pca,
    run_deseq2,
    significant,
)


def synthetic(seed: int = 0, genes: int = 400) -> RnaSeqDataset:
    rng = np.random.default_rng(seed)
    mu = rng.gamma(2.0, 200.0, size=genes)
    depth = np.array([0.5, 1.0, 2.0, 0.8, 1.2, 1.6])
    effect = np.ones(genes)
    effect[:40] = 8.0  # 40 genes up 8x in group B
    group = np.array([0, 0, 0, 1, 1, 1])
    lam = mu[:, None] * depth[None, :] * np.where(group[None, :] == 1, effect[:, None], 1.0)
    counts = rng.negative_binomial(n=20, p=20 / (20 + lam))
    samples = pd.DataFrame({"condition": ["A"] * 3 + ["B"] * 3}, index=[f"s{i}" for i in range(6)])
    return RnaSeqDataset(
        "toy",
        pd.DataFrame(counts, index=[f"g{i}" for i in range(genes)], columns=samples.index),
        samples,
    )


def test_size_factors_recover_known_depths():
    base = np.arange(1, 101, dtype=float)[:, None]
    depth = np.array([1.0, 2.0, 4.0])
    counts = pd.DataFrame(base * depth, columns=["a", "b", "c"])
    sf = median_of_ratios_size_factors(counts)
    np.testing.assert_allclose(sf / sf.iloc[0], depth)
    assert np.exp(np.log(sf).mean()) == pytest.approx(1.0)
    np.testing.assert_allclose(normalized_counts(counts).std(axis=1), 0, atol=1e-9)


def test_size_factors_need_a_gene_expressed_everywhere():
    with pytest.raises(ValueError, match="undefined"):
        median_of_ratios_size_factors(pd.DataFrame({"a": [0, 1], "b": [1, 0]}))


def test_dataset_rejects_mismatched_samples():
    with pytest.raises(ValueError, match="disagree"):
        RnaSeqDataset(
            "x", pd.DataFrame({"a": [1]}), pd.DataFrame({"condition": ["A"]}, index=["b"])
        )


def test_pca_separates_groups_and_is_sign_stable():
    ds = synthetic()
    coords, explained = pca(log_normalized(ds.counts), n_top=200)
    assert explained[0] > 0.5
    assert np.sign(coords.loc[["s0", "s1", "s2"], "PC1"]).nunique() == 1
    assert np.sign(coords.loc["s0", "PC1"]) != np.sign(coords.loc["s3", "PC1"])


def test_deseq2_on_synthetic_data_finds_planted_genes():
    ds = synthetic()
    res, sf = run_deseq2(ds, "condition", reference="A", test="B")
    np.testing.assert_allclose(sf, median_of_ratios_size_factors(ds.counts), rtol=1e-10)
    hits = set(significant(res, 0.05, 1.0).index)
    planted = {f"g{i}" for i in range(40)}
    assert len(hits & planted) >= 36 and len(hits - planted) <= 2
    assert direction_counts(res)["up"] >= 36
    with pytest.raises(ValueError, match="not in design"):
        run_deseq2(ds, "condition", reference="A", test="Z")


def test_compare_results_identity():
    ds = synthetic(1)
    res, _ = run_deseq2(ds, "condition", "A", "B")
    m = compare_results(res, res.copy())
    assert (
        m["lfc_pearson"] == pytest.approx(1.0)
        and m["jaccard"] == 1.0
        and m["lfc_max_abs_diff"] == 0
    )


def test_normalized_counts_equal_galaxy_deseq2(data_dir):
    ds = load_dataset("mabs", data_dir / "rnaseq/mabs/counts", data_dir / "rnaseq/mabs/samples.tsv")
    galaxy = pd.read_csv(
        data_dir / "rnaseq/mabs/galaxy_normalized_counts.tsv", sep="\t", index_col=0
    )
    ours = normalized_counts(ds.counts)[galaxy.columns].loc[galaxy.index]
    np.testing.assert_allclose(ours.to_numpy(), galaxy.to_numpy(), rtol=1e-12)


@pytest.mark.slow
@pytest.mark.parametrize(
    ("name", "reference", "test", "min_recovered"),
    [("mabs", "WT_rough", "condition", 0.99), ("mtb", "control", "experimental", 1.0)],
)
def test_pydeseq2_reproduces_galaxy_deseq2(data_dir, name, reference, test, min_recovered):
    ds = load_dataset(
        name, data_dir / f"rnaseq/{name}/counts", data_dir / f"rnaseq/{name}/samples.tsv"
    )
    ours, _ = run_deseq2(ds, "condition", reference, test)
    galaxy = read_deseq2_table(data_dir / f"rnaseq/{name}/galaxy_deseq2_results.tsv")
    m = compare_results(ours, galaxy)
    assert m["basemean_pearson"] > 0.999999
    assert m["neglog10p_spearman"] > 0.999
    assert m["significant_both"] / m["significant_reference"] >= min_recovered
    assert m["sign_agreement_pct"] == 100.0
