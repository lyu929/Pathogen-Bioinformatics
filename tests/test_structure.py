import numpy as np
import pytest

from pathogenomics.structure import (
    ChainCA,
    af3_pae,
    af3_query_sequence,
    af3_summary,
    align_sequences,
    chain_ca,
    compare_chains,
    iterative_core_rmsd,
    kabsch,
    load_structure,
    plddt_bands,
    read_seqres,
    residue_correspondence,
    rigid_domains,
)


def random_rotation(rng) -> np.ndarray:
    qm, r = np.linalg.qr(rng.normal(size=(3, 3)))
    qm = qm @ np.diag(np.sign(np.diag(r)))
    if np.linalg.det(qm) < 0:
        qm[:, 0] *= -1
    return qm


def test_kabsch_recovers_rigid_motion():
    rng = np.random.default_rng(0)
    x = rng.normal(size=(50, 3)) * 10
    rot = random_rotation(rng)
    y = x @ rot.T + np.array([3.0, -2.0, 7.0])
    r, t, rmsd = kabsch(x, y)
    assert rmsd == pytest.approx(0, abs=1e-9)
    np.testing.assert_allclose(r, rot, atol=1e-9)
    assert np.linalg.det(r) == pytest.approx(1.0)


def test_kabsch_never_returns_a_reflection():
    rng = np.random.default_rng(1)
    x = rng.normal(size=(30, 3))
    mirrored = x * np.array([1, 1, -1])
    r, _, rmsd = kabsch(x, mirrored)
    assert np.linalg.det(r) == pytest.approx(1.0) and rmsd > 0.1
    with pytest.raises(ValueError):
        kabsch(x[:2], x[:2])


def test_core_rmsd_ignores_a_floppy_tail():
    rng = np.random.default_rng(2)
    x = rng.normal(size=(60, 3)) * 8
    y = x.copy()
    y[50:] += 15.0
    rmsd, n = iterative_core_rmsd(x, y)
    assert n == 50 and rmsd == pytest.approx(0, abs=1e-6)


def make_chain(coords: np.ndarray, seq: str, start: int = 1) -> ChainCA:
    return ChainCA(np.arange(start, start + len(seq)), seq, coords, np.full(len(seq), 90.0))


def test_rigid_domains_detects_a_hinge():
    rng = np.random.default_rng(3)
    d1, d2 = rng.normal(size=(40, 3)) * 5, rng.normal(size=(40, 3)) * 5 + 30
    a = np.vstack([d1, d2])
    rot = random_rotation(rng)
    b = np.vstack([d1, (d2 - d2.mean(0)) @ rot.T + d2.mean(0) + 5])  # domain 2 swings
    seq = "A" * 80
    ca, cb = make_chain(a, seq), make_chain(b, seq)
    pairs, _ = residue_correspondence(ca, cb)
    bodies = rigid_domains(ca, cb, pairs)
    assert len(bodies) == 2 and (bodies["pairs"] == 40).all() and (bodies["rmsd_ca"] < 1e-6).all()
    assert compare_chains(ca, cb)["rmsd_all_ca"] > 3


def test_alignment_identity_and_seqres_projection():
    aln = align_sequences("MKTAYIAKQR", "MKTAYLAKQR")
    assert len(aln.pairs) == 10 and aln.identity == pytest.approx(0.9)
    # chain b is missing residues 3-4 (unmodelled) and has one substitution at 8
    seqres_a, seqres_b = "MKTAYIAKQRGG", "MKTAYIAWQRGG"
    rng = np.random.default_rng(4)
    xyz = rng.normal(size=(12, 3))
    a = make_chain(xyz, seqres_a)
    keep = [i for i in range(12) if i not in (2, 3)]
    b = ChainCA(
        np.arange(1, 13)[keep], "".join(seqres_b[i] for i in keep), xyz[keep], np.full(10, 90.0)
    )
    pairs, diffs = residue_correspondence(a, b, seqres_a, seqres_b)
    assert len(pairs) == 10 and diffs["label"].tolist() == ["K8W"]


def test_plddt_bands():
    bands = plddt_bands(np.array([95, 80, 60, 10]))
    assert bands == {
        "very_high": 0.25,
        "confident": 0.25,
        "low": 0.25,
        "very_low": 0.25,
        "mean_plddt": 61.25,
    }


def test_course_alphafold3_zur(data_dir):
    d = data_dir / "structure/alphafold3_zur_Rv2359"
    summary = af3_summary(d)
    assert summary["ranking_score"].idxmax() == 0 and summary.loc[0, "ptm"] == pytest.approx(0.75)
    model = chain_ca(load_structure(d / "model_0.cif"))
    assert model.sequence == af3_query_sequence(d / "job_request.json")
    assert af3_pae(d / "full_data_0.json").shape == (len(model.sequence),) * 2


def test_course_spike_rbd_mutations(data_dir):
    a = chain_ca(load_structure(data_dir / "structure/pdb/6vxx.pdb.gz"), "A")
    b = chain_ca(load_structure(data_dir / "structure/pdb/9f9y.pdb.gz"), "A")
    _, diffs = residue_correspondence(
        a,
        b,
        read_seqres(data_dir / "structure/pdb/6vxx.pdb.gz")["A"],
        read_seqres(data_dir / "structure/pdb/9f9y.pdb.gz")["A"],
    )
    labels = set(diffs["label"])
    assert {"D614G", "N501Y", "K417T", "S373P", "N440K"} <= labels


def test_course_dosr_homologs_are_two_rigid_domains(data_dir):
    a = chain_ca(load_structure(data_dir / "structure/afdb/AF-B1MH84-F1-model_v6.pdb"))
    b = chain_ca(load_structure(data_dir / "structure/afdb/AF-P9WMF9-F1-model_v6.pdb"))
    pairs, _ = residue_correspondence(a, b)
    bodies = rigid_domains(a, b, pairs)
    assert len(bodies) == 2 and (bodies["rmsd_ca"] < 1.0).all()
    assert compare_chains(a, b)["identity_pct"] == pytest.approx(70.4, abs=0.1)
