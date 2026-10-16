"""Protein structure analysis: AlphaFold confidence metrics and sequence-guided superposition."""

from __future__ import annotations

import gzip
import io
import json
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}  # fmt: skip

PLDDT_BANDS = (("very_high", 90.0), ("confident", 70.0), ("low", 50.0), ("very_low", 0.0))


# --------------------------------------------------------------------------- AlphaFold 3
def af3_summary(directory: str | Path) -> pd.DataFrame:
    """Per-model AlphaFold Server summary confidences (``summary_confidences_*.json``)."""
    rows = []
    for path in sorted(Path(directory).glob("summary_confidences_*.json")):
        data = json.loads(path.read_text())
        rows.append(
            {
                "model": int(path.stem.rsplit("_", 1)[-1]),
                "ranking_score": data["ranking_score"],
                "ptm": data["ptm"],
                "iptm": data.get("iptm"),
                "fraction_disordered": data["fraction_disordered"],
                "has_clash": bool(data["has_clash"]),
                "num_recycles": int(data.get("num_recycles", 0)),
            }
        )
    if not rows:
        raise FileNotFoundError(f"no summary_confidences_*.json in {directory}")
    return pd.DataFrame(rows).set_index("model").sort_index()


def af3_pae(full_data_json: str | Path) -> np.ndarray:
    """Predicted aligned error matrix (tokens x tokens, Angstrom) from ``full_data_*.json``."""
    data = json.loads(Path(full_data_json).read_text())
    pae = np.asarray(data["pae"], dtype=float)
    if pae.ndim != 2 or pae.shape[0] != pae.shape[1]:
        raise ValueError("PAE matrix must be square")
    return pae


def af3_query_sequence(job_request_json: str | Path) -> str:
    job = json.loads(Path(job_request_json).read_text())
    job = job[0] if isinstance(job, list) else job
    chains = [s["proteinChain"]["sequence"] for s in job["sequences"] if "proteinChain" in s]
    if len(chains) != 1:
        raise ValueError("expected a single protein chain in the AlphaFold job request")
    return chains[0]


# --------------------------------------------------------------------------- structure parsing
def _open_structure_text(path: Path) -> io.StringIO:
    raw = path.read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return io.StringIO(raw.decode("utf-8"))


def load_structure(path: str | Path, structure_id: str | None = None):
    """Parse a PDB or mmCIF file (optionally gzipped) with Biopython."""
    from Bio.PDB import MMCIFParser, PDBParser

    path = Path(path)
    name = path.name.removesuffix(".gz")
    parser = MMCIFParser(QUIET=True) if name.endswith(".cif") else PDBParser(QUIET=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return parser.get_structure(structure_id or name.split(".")[0], _open_structure_text(path))


@dataclass
class ChainCA:
    """CA trace of one chain: residue numbers, one-letter sequence, coordinates, B-factors."""

    resnums: np.ndarray
    sequence: str
    coords: np.ndarray
    bfactors: np.ndarray


def chain_ca(structure, chain_id: str | None = None, model: int = 0) -> ChainCA:
    model_obj = list(structure)[model]
    chain = model_obj[chain_id] if chain_id is not None else next(iter(model_obj))
    nums, seq, xyz, bf = [], [], [], []
    for residue in chain:
        hetflag, resnum, _icode = residue.id
        if hetflag != " " or "CA" not in residue:
            continue
        aa = THREE_TO_ONE.get(residue.get_resname())
        if aa is None:
            continue
        atom = residue["CA"]
        nums.append(resnum)
        seq.append(aa)
        xyz.append(atom.coord)
        bf.append(atom.bfactor)
    if not seq:
        raise ValueError(f"chain {chain.id} has no standard amino-acid CA atoms")
    return ChainCA(
        np.asarray(nums), "".join(seq), np.asarray(xyz, dtype=float), np.asarray(bf, dtype=float)
    )


def plddt_bands(plddt: np.ndarray) -> dict[str, float]:
    """Fraction of residues in AlphaFold's pLDDT confidence bands."""
    plddt = np.asarray(plddt, dtype=float)
    out: dict[str, float] = {}
    upper = np.inf
    for name, lower in PLDDT_BANDS:
        out[name] = round(float(((plddt >= lower) & (plddt < upper)).mean()), 4)
        upper = lower
    out["mean_plddt"] = round(float(plddt.mean()), 2)
    return out


# --------------------------------------------------------------------------- alignment + superposition
@dataclass
class SequenceAlignment:
    pairs: np.ndarray  # (k, 2) indices into seq_a / seq_b of aligned (non-gap) columns
    identity: float  # identical / aligned columns
    coverage_a: float
    coverage_b: float
    score: float


def align_sequences(seq_a: str, seq_b: str) -> SequenceAlignment:
    """Global protein alignment (BLOSUM62, gap open -10, extend -0.5) with free end gaps."""
    from Bio import Align
    from Bio.Align import substitution_matrices

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    aligner.end_gap_score = 0.0
    aligner.mode = "global"
    best = aligner.align(seq_a, seq_b)[0]
    pairs = []
    for (a0, a1), (b0, b1) in zip(*best.aligned, strict=True):
        pairs.extend(zip(range(a0, a1), range(b0, b1), strict=True))
    pairs_arr = np.asarray(pairs, dtype=int).reshape(-1, 2)
    identical = sum(seq_a[i] == seq_b[j] for i, j in pairs_arr)
    n = max(len(pairs_arr), 1)
    return SequenceAlignment(
        pairs=pairs_arr,
        identity=identical / n,
        coverage_a=len(pairs_arr) / len(seq_a),
        coverage_b=len(pairs_arr) / len(seq_b),
        score=float(best.score),
    )


def kabsch(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Optimal rotation ``R`` and translation ``t`` minimising RMSD of ``mobile @ R.T + t`` to ``target``.

    Returns ``(R, t, rmsd)``. Reflections are excluded (proper rotation, det(R) = +1).
    """
    mobile = np.asarray(mobile, dtype=float)
    target = np.asarray(target, dtype=float)
    if mobile.shape != target.shape or mobile.ndim != 2 or mobile.shape[1] != 3:
        raise ValueError("coordinate arrays must both be (n, 3)")
    if mobile.shape[0] < 3:
        raise ValueError("need at least 3 points to superpose")
    mu_m, mu_t = mobile.mean(axis=0), target.mean(axis=0)
    h = (mobile - mu_m).T @ (target - mu_t)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    trans = mu_t - mu_m @ rot.T
    moved = mobile @ rot.T + trans
    rmsd = float(np.sqrt(((moved - target) ** 2).sum(axis=1).mean()))
    return rot, trans, rmsd


def _refine(
    mobile: np.ndarray, target: np.ndarray, keep: np.ndarray, cutoff: float, max_iter: int = 30
) -> np.ndarray:
    """Re-fit on the current inliers and re-select pairs within ``cutoff`` until stable."""
    for _ in range(max_iter):
        if keep.sum() < 3:
            return keep
        rot, trans, _ = kabsch(mobile[keep], target[keep])
        new_keep = np.linalg.norm(mobile @ rot.T + trans - target, axis=1) <= cutoff
        if np.array_equal(new_keep, keep):
            break
        keep = new_keep
    return keep


def largest_rigid_subset(
    mobile: np.ndarray, target: np.ndarray, cutoff: float = 2.0, window: int = 12, n_refine: int = 5
) -> np.ndarray:
    """Boolean mask of the largest set of pairs that superpose within ``cutoff`` Angstrom.

    Every contiguous ``window`` of pairs seeds a superposition; the seeds with the most
    inliers are refined (fit on inliers -> re-select) and the largest consensus set wins.
    Seeding from local fragments, rather than from a global fit, keeps the search from
    being dragged off by a second domain or a displaced tail.
    """
    n = len(mobile)
    if n < 3:
        return np.zeros(n, dtype=bool)
    window = max(3, min(window, n))
    scores = []
    for start in range(0, n - window + 1):
        sl = slice(start, start + window)
        rot, trans, _ = kabsch(mobile[sl], target[sl])
        inliers = np.linalg.norm(mobile @ rot.T + trans - target, axis=1) <= cutoff
        scores.append((int(inliers.sum()), start, inliers))
    scores.sort(key=lambda x: (-x[0], x[1]))
    best = np.zeros(n, dtype=bool)
    for _, _, inliers in scores[:n_refine]:
        refined = _refine(mobile, target, inliers, cutoff)
        if refined.sum() > best.sum():
            best = refined
    return best


def iterative_core_rmsd(
    mobile: np.ndarray, target: np.ndarray, cutoff: float = 2.0
) -> tuple[float, int]:
    """RMSD over the largest rigid core (pairs within ``cutoff`` after superposition).

    Returns ``(core_rmsd, n_core_pairs)`` - a measure robust to flexible loops, termini
    and hinge motions, in the spirit of PyMOL ``align`` outlier rejection.
    """
    keep = largest_rigid_subset(mobile, target, cutoff)
    if keep.sum() < 3:
        return float("nan"), int(keep.sum())
    return kabsch(mobile[keep], target[keep])[2], int(keep.sum())


def read_seqres(path: str | Path) -> dict[str, str]:
    """Full deposited sequence per chain from PDB ``SEQRES`` records (includes unmodelled residues)."""
    from Bio import SeqIO

    handle = _open_structure_text(Path(path))
    return {rec.id.split(":")[-1]: str(rec.seq) for rec in SeqIO.parse(handle, "pdb-seqres")}


def _project_to_seqres(chain: ChainCA, seqres: str) -> np.ndarray:
    """Index into ``seqres`` of every modelled residue (-1 if it cannot be placed)."""
    out = np.full(len(chain.sequence), -1, dtype=int)
    aln = align_sequences(chain.sequence, seqres)
    for i, j in aln.pairs:
        if chain.sequence[i] == seqres[j]:
            out[i] = j
    return out


def residue_correspondence(
    a: ChainCA, b: ChainCA, seqres_a: str | None = None, seqres_b: str | None = None
) -> tuple[np.ndarray, pd.DataFrame]:
    """Pairs of corresponding modelled residues, plus differences in ``a``'s numbering.

    Aligning only the *modelled* residues of two chains misplaces residues next to
    unresolved loops. When the full ``SEQRES`` sequences are given, they are aligned
    instead and each modelled residue is projected onto its own SEQRES, so
    correspondences (and substitutions/deletions) are read off the full-length alignment.

    Returns ``(pairs, differences)`` where ``pairs`` is a (k, 2) array of indices into
    ``a``/``b`` and ``differences`` lists substitutions and deletions at residues
    modelled in ``a``.
    """
    rows: list[dict[str, object]] = []
    if seqres_a is None or seqres_b is None:
        aln = align_sequences(a.sequence, b.sequence)
        pairs = aln.pairs
        for i, j in pairs:
            if a.sequence[i] != b.sequence[j]:
                rows.append(
                    {"position": int(a.resnums[i]), "from": a.sequence[i], "to": b.sequence[j]}
                )
    else:
        pa, pb = _project_to_seqres(a, seqres_a), _project_to_seqres(b, seqres_b)
        full = align_sequences(seqres_a, seqres_b)
        sa_to_sb = dict(zip(full.pairs[:, 0], full.pairs[:, 1], strict=True))
        sb_to_b = {int(s): k for k, s in enumerate(pb) if s >= 0}
        pair_list = []
        first, last = full.pairs[0, 0], full.pairs[-1, 0]
        for i, s_a in enumerate(pa):
            if s_a < 0:
                continue
            s_b = sa_to_sb.get(int(s_a))
            if s_b is None:
                if first <= s_a <= last:  # internal gap in b -> deletion
                    rows.append({"position": int(a.resnums[i]), "from": a.sequence[i], "to": "-"})
                continue
            if seqres_b[s_b] != a.sequence[i]:
                rows.append(
                    {"position": int(a.resnums[i]), "from": a.sequence[i], "to": seqres_b[s_b]}
                )
            if s_b in sb_to_b:
                pair_list.append((i, sb_to_b[s_b]))
        pairs = np.asarray(pair_list, dtype=int).reshape(-1, 2)
    diffs = pd.DataFrame(rows, columns=["position", "from", "to"])
    diffs["label"] = [
        f"{f}{p}{t}" if t != "-" else f"del{f}{p}"
        for p, f, t in zip(diffs["position"], diffs["from"], diffs["to"], strict=True)
    ]
    return pairs, diffs


def compare_chains(
    a: ChainCA,
    b: ChainCA,
    seqres_a: str | None = None,
    seqres_b: str | None = None,
    plddt_min: float | None = None,
) -> dict[str, float | int]:
    """Superpose corresponding CA atoms of two chains (see :func:`residue_correspondence`).

    ``plddt_min`` (for AlphaFold models, whose B-factor column holds pLDDT) adds an RMSD
    restricted to residue pairs where both models are at least that confident.
    """
    pairs, diffs = residue_correspondence(a, b, seqres_a, seqres_b)
    ia, ib = pairs[:, 0], pairs[:, 1]
    identical = int(sum(a.sequence[i] == b.sequence[j] for i, j in pairs))
    _, _, rmsd_all = kabsch(a.coords[ia], b.coords[ib])
    core_rmsd, core_n = iterative_core_rmsd(a.coords[ia], b.coords[ib])
    out: dict[str, float | int] = {
        "modelled_a": len(a.sequence),
        "modelled_b": len(b.sequence),
        "aligned_pairs": len(pairs),
        "identity_pct": round(100 * identical / max(len(pairs), 1), 1),
        "substitutions": int((diffs["to"] != "-").sum()),
        "deletions": int((diffs["to"] == "-").sum()),
        "rmsd_all_ca": round(rmsd_all, 2),
        "rmsd_core_ca": round(core_rmsd, 2),
        "core_pairs": core_n,
    }
    if plddt_min is not None:
        conf = (a.bfactors[ia] >= plddt_min) & (b.bfactors[ib] >= plddt_min)
        if conf.sum() >= 3:
            out[f"rmsd_ca_plddt{int(plddt_min)}"] = round(
                kabsch(a.coords[ia][conf], b.coords[ib][conf])[2], 2
            )
            out[f"pairs_plddt{int(plddt_min)}"] = int(conf.sum())
    return out


def rigid_domains(
    a: ChainCA,
    b: ChainCA,
    pairs: np.ndarray,
    cutoff: float = 2.0,
    min_size: int = 20,
    max_domains: int = 3,
) -> pd.DataFrame:
    """Decompose a superposition into rigid bodies by repeated core finding.

    The largest set of residue pairs that superposes within ``cutoff`` Angstrom
    (:func:`largest_rigid_subset`) is taken as the first rigid body, removed, and the
    search repeated on the remaining pairs. Low per-body RMSD together with a high global
    RMSD means the chains share domain structures but differ in domain orientation.
    """
    remaining = np.ones(len(pairs), dtype=bool)
    rows = []
    for k in range(max_domains):
        idx = np.flatnonzero(remaining)
        if idx.size < min_size:
            break
        keep = largest_rigid_subset(a.coords[pairs[idx, 0]], b.coords[pairs[idx, 1]], cutoff)
        if keep.sum() < min_size:
            break
        members = idx[keep]
        rmsd = kabsch(a.coords[pairs[members, 0]], b.coords[pairs[members, 1]])[2]
        resnums = a.resnums[pairs[members, 0]]
        rows.append(
            {
                "body": k + 1,
                "pairs": int(members.size),
                "rmsd_ca": round(rmsd, 2),
                "first_residue": int(resnums.min()),
                "last_residue": int(resnums.max()),
                "median_residue": int(np.median(resnums)),
            }
        )
        remaining[members] = False
    return pd.DataFrame(
        rows,
        columns=["body", "pairs", "rmsd_ca", "first_residue", "last_residue", "median_residue"],
    )
