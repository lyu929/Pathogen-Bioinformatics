import numpy as np
import pandas as pd
import pytest

from pathogenomics.pangenome import (
    ROARY_METADATA,
    absent_in,
    infer_isolate_names,
    jaccard_distance,
    partition,
    partition_summary,
    presence_matrix,
    rarefaction,
    read_roary,
    unique_genes,
)


def toy_roary(tmp_path):
    rows = [
        ("core1", ["A_00001", "B_00001", "C_00001", "D_00001"]),
        ("core2", ["A_00002", "B_00002", "C_00002", "D_00002"]),
        ("shell1", ["A_00003", "B_00003", "", ""]),
        ("cloudA", ["A_00004", "", "", ""]),
        ("missD", ["A_00005", "B_00005", "C_00005", ""]),
    ]
    records = []
    for gene, cells in rows:
        meta = {c: "" for c in ROARY_METADATA}
        meta["Gene"] = gene
        meta["Annotation"] = f"{gene} protein"
        meta["No. isolates"] = str(sum(bool(c) for c in cells))
        records.append(
            {
                **meta,
                **dict(zip(["Prokka_1", "Prokka_2", "Prokka_3", "Prokka_4"], cells, strict=True)),
            }
        )
    path = tmp_path / "gpa.csv"
    pd.DataFrame(records).to_csv(path, index=False)
    return read_roary(path)


def test_isolate_names_from_locus_prefixes(tmp_path):
    table = toy_roary(tmp_path)
    assert infer_isolate_names(table) == {
        "Prokka_1": "A",
        "Prokka_2": "B",
        "Prokka_3": "C",
        "Prokka_4": "D",
    }


def test_partition_boundaries(tmp_path):
    p = presence_matrix(toy_roary(tmp_path))
    labels = partition(p)
    assert labels.to_dict() == {
        "core1": "core",
        "core2": "core",
        "shell1": "shell",
        "cloudA": "shell",
        "missD": "shell",
    }
    # with 20 genomes a singleton (5%) is cloud and 19/20 (95%) is soft core
    big = pd.DataFrame(np.ones((3, 20), dtype=bool), index=["c", "s", "x"])
    big.loc["s", big.columns[:1]] = False
    big.loc["x"] = False
    big.loc["x", 0] = True
    assert partition(big).to_dict() == {"c": "core", "s": "soft_core", "x": "cloud"}
    summary = partition_summary(p)
    assert summary.loc["total", "clusters"] == 5 and summary["clusters"].dtype.kind == "i"


def test_rarefaction_monotone_and_deterministic(tmp_path):
    p = presence_matrix(toy_roary(tmp_path))
    a = rarefaction(p, permutations=50, seed=3)
    b = rarefaction(p, permutations=50, seed=3)
    pd.testing.assert_frame_equal(a, b)
    assert (np.diff(a["pan_mean"]) >= 0).all() and (np.diff(a["core_mean"]) <= 0).all()
    assert a["pan_mean"].iloc[-1] == 5 and a["core_mean"].iloc[-1] == 2


def test_distances_and_unique_genes(tmp_path):
    table = toy_roary(tmp_path)
    p = presence_matrix(table)
    d = jaccard_distance(p)
    assert np.allclose(d.values, d.values.T) and (np.diag(d.values) == 0).all()
    assert d.loc["A", "B"] == pytest.approx(1 - 2 / 3)  # accessory: shell1, cloudA, missD
    assert unique_genes(table, p).to_dict("index") == {
        "cloudA": {"isolate": "A", "annotation": "cloudA protein"}
    }
    assert list(absent_in(p, table, "D").index) == ["missD"]


def test_presence_must_match_roary_counts(tmp_path):
    table = toy_roary(tmp_path)
    table.loc["core1", "No. isolates"] = 3
    with pytest.raises(ValueError, match="disagrees"):
        presence_matrix(table)


def test_course_roary_output(data_dir):
    table = read_roary(data_dir / "pangenome/roary_gene_presence_absence.csv.gz")
    p = presence_matrix(table)
    assert list(p.columns) == ["EN93", "EN94", "EN96", "EN97", "PA14"]
    s = partition_summary(p)
    assert s.loc["core", "clusters"] == 5874 and s.loc["total", "clusters"] == 5938
