"""End-to-end pipeline: every stage writes tables + figures and returns metrics for the report."""

from __future__ import annotations

import dataclasses
import json
import logging
import time
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from pathogenomics import annotation, assembly, pangenome, plotting, qc, rnaseq, structure
from pathogenomics.io.galaxy import read_deseq2_table
from pathogenomics.io.gff import gene_annotation_table
from pathogenomics.io.reports import read_cutadapt_report, read_fastqc_data

log = logging.getLogger("pathogenomics")

STAGES = ("qc", "assembly", "annotation", "pangenome", "rnaseq", "structure")


class Context:
    def __init__(self, cfg: dict[str, Any]):
        self.cfg = cfg
        self.tables = Path(cfg["output"]["results"])
        self.report_dir = Path(cfg["output"]["report"])
        self.figures = self.report_dir / "figures"
        self.tables.mkdir(parents=True, exist_ok=True)
        self.figures.mkdir(parents=True, exist_ok=True)

    def table(self, df: pd.DataFrame, name: str, index: bool = True) -> Path:
        path = self.tables / f"{name}.tsv"
        df.to_csv(path, sep="\t", index=index, float_format="%.6g")
        return path

    def figure(self, draw: Callable, stem: str) -> str:
        plotting.save_both(draw, self.figures, stem)
        return stem


# --------------------------------------------------------------------------- stages
def run_qc(ctx: Context) -> dict[str, Any]:
    c = ctx.cfg["qc"]
    trimmed_path = ctx.tables / "BR-Tn1.trimmed.fastq.gz"
    summary = qc.trim_fastq(c["fastq"], c["quality_cutoff"], c["min_length"], output=trimmed_path)
    galaxy = read_cutadapt_report(c["cutadapt_of_fastq"])
    raw_stats = qc.fastq_stats(c["fastq"])
    trimmed_stats = qc.fastq_stats(trimmed_path)
    trimmed_path.unlink()  # derived data; regenerated on every run

    fastqc_ref = read_fastqc_data(c["fastqc_of_fastq"])
    bins = fastqc_ref.per_base_quality()
    ours_binned = qc.binned_mean_quality(raw_stats, bins)
    max_q_diff = float(np.abs(ours_binned - bins["Mean"].to_numpy()).max())

    fastqc_reports = {
        p.name.split(".")[0]: read_fastqc_data(p)
        for p in sorted(Path(c["fastqc_dir"]).glob("*.txt*"))
    }
    cutadapt_reports = {
        p.name.split(".")[0]: read_cutadapt_report(p)
        for p in sorted(Path(c["cutadapt_dir"]).glob("*.txt"))
    }
    fastqc_table = qc.fastqc_summary(fastqc_reports)
    fastqc_table["failed_modules"] = qc.flag_failing_modules(fastqc_table)
    cut_table = qc.cutadapt_summary(cutadapt_reports)
    ctx.table(fastqc_table, "qc_fastqc_summary")
    ctx.table(cut_table, "qc_cutadapt_summary")

    ours = dataclasses.asdict(summary)
    validation = pd.DataFrame(
        {"galaxy_cutadapt": [getattr(galaxy, k) for k in ours], "this_repo": list(ours.values())},
        index=pd.Index(list(ours), name="count"),
    )
    validation["match"] = validation["galaxy_cutadapt"] == validation["this_repo"]
    ctx.table(validation, "qc_trimming_validation")

    fig = ctx.figure(
        plotting.per_base_quality(
            {
                "raw reads": raw_stats.per_position_mean_quality,
                "after Q20 trimming": trimmed_stats.per_position_mean_quality,
            },
            {
                "raw reads": raw_stats.per_position_count,
                "after Q20 trimming": trimmed_stats.per_position_count,
            },
            c["quality_cutoff"],
        ),
        "qc_per_base_quality",
    )
    return {
        "trim_matches_cutadapt": summary.matches(galaxy),
        "trim_validation": validation,
        "cutadapt_version": galaxy.version,
        "raw": raw_stats.summary(),
        "trimmed": trimmed_stats.summary(),
        "fastqc_reads": fastqc_ref.total_sequences,
        "fastqc_gc": fastqc_ref.gc_percent,
        "fastqc_max_mean_quality_diff": round(max_q_diff, 6),
        "fastqc_table": fastqc_table,
        "cutadapt_table": cut_table,
        "figure": fig,
    }


def run_assembly(ctx: Context) -> dict[str, Any]:
    c = ctx.cfg["assembly"]
    isolates = sorted(Path(c["isolates_dir"]).glob("*.fasta.gz"))
    table = assembly.assembly_table(isolates, min_length=c["min_contig"])
    table.index = [i.split(".")[0] for i in table.index]
    refs = pd.read_csv(c["reference_genomes"], sep="\t")
    other = assembly.assembly_table(c["other"], min_length=0)
    ctx.table(table, "assembly_pseudomonas_isolates")
    ctx.table(other, "assembly_other")
    lengths = {p.name.split(".")[0]: assembly.contig_lengths(p, c["min_contig"]) for p in isolates}
    fig = ctx.figure(plotting.nx_curves(lengths), "assembly_nx_curves")
    return {
        "table": table,
        "other": other,
        "ref_min": int(refs["size_bp"].min()),
        "ref_max": int(refs["size_bp"].max()),
        "ref_n": len(refs),
        "figure": fig,
    }


def run_annotation(ctx: Context) -> dict[str, Any]:
    c = ctx.cfg["annotation"]
    strains = pd.read_csv(c["strains"], sep="\t").set_index("accession")
    gffs: dict[str, Path] = {}
    for path in sorted(Path(c["refseq_dir"]).glob("*.gff.gz")):
        acc = path.name.split(".gff")[0]
        gffs[strains.loc[acc, "strain"]] = path
    rows = [annotation.annotation_summary(p, name) for name, p in gffs.items()]
    genome_len = assembly.assembly_stats(c["bvbrc_genome"]).total_length
    rows.append(
        annotation.annotation_summary(c["bvbrc_gff"], c["bvbrc_label"], genome_length=genome_len)
    )
    summary = pd.DataFrame(rows).set_index("name")
    gffs[c["bvbrc_label"]] = c["bvbrc_gff"]
    presence, evidence = annotation.marker_matrix(gffs)
    collisions = annotation.symbol_collisions(gffs)
    collision_summary = (
        collisions.groupby(["symbol", "product"]).genome.nunique().rename("genomes").reset_index()
        if not collisions.empty
        else collisions
    )
    ctx.table(summary, "annotation_summary")
    ctx.table(presence.astype(int), "annotation_marker_presence")
    ctx.table(evidence, "annotation_marker_evidence", index=False)
    ctx.table(collision_summary, "annotation_symbol_collisions", index=False)
    fig = ctx.figure(
        plotting.presence_heatmap(
            presence,
            "Resistance and virulence markers (from protein product names)",
            cell_labels=True,
        ),
        "annotation_markers",
    )
    mrsa = presence.index[presence["mecA (PBP2a)"]].tolist()
    return {
        "summary": summary,
        "presence": presence,
        "collisions": collision_summary,
        "mrsa": mrsa,
        "mssa": [g for g in presence.index if g not in mrsa],
        "figure": fig,
    }


def run_pangenome(ctx: Context) -> dict[str, Any]:
    c = ctx.cfg["pangenome"]
    table = pangenome.read_roary(c["roary"])
    presence = pangenome.presence_matrix(table)
    parts = pangenome.partition_summary(presence)
    curve = pangenome.rarefaction(presence, permutations=c["permutations"], seed=c["seed"])
    dist = pangenome.jaccard_distance(presence)
    uniq = pangenome.unique_genes(table, presence)
    accessory = presence[~presence.all(axis=1)]
    ctx.table(parts, "pangenome_partitions")
    ctx.table(curve, "pangenome_rarefaction", index=False)
    ctx.table(dist.round(4), "pangenome_jaccard_accessory")
    ctx.table(uniq, "pangenome_unique_genes")
    ordered = accessory.loc[accessory.sum(axis=1).sort_values(ascending=False).index]
    fig_curve = ctx.figure(plotting.rarefaction_plot(curve), "pangenome_rarefaction")
    fig_heat = ctx.figure(
        plotting.presence_heatmap(
            ordered.T,
            f"Accessory gene clusters ({len(accessory)}) across isolates",
            xlabel="Roary gene cluster",
        ),
        "pangenome_accessory",
    )
    isolates = [c for c in presence.columns if c != "PA14"]
    iso_dist = dist.loc[isolates, isolates].to_numpy()
    return {
        "partitions": parts,
        "n_clusters": len(presence),
        "n_genomes": presence.shape[1],
        "accessory": len(accessory),
        "unique_counts": uniq["isolate"].value_counts(),
        "distance": dist,
        "isolate_dist_range": (
            float(iso_dist[np.triu_indices(len(isolates), 1)].min()),
            float(iso_dist[np.triu_indices(len(isolates), 1)].max()),
        ),
        "reference_dist_min": float(dist.loc[isolates, "PA14"].min()) if "PA14" in dist else None,
        "figures": [fig_curve, fig_heat],
    }


def run_rnaseq(ctx: Context) -> dict[str, Any]:
    th = ctx.cfg["thresholds"]
    out: dict[str, Any] = {}
    for c in ctx.cfg["rnaseq"]:
        name = c["name"]
        ds = rnaseq.load_dataset(name, c["counts_dir"], c["samples"])
        ann = gene_annotation_table(c["annotation"])
        lib = rnaseq.library_summary(ds)
        biotype = ann["biotype"].reindex(ds.counts.index).fillna("unknown")
        lib["rRNA_pct"] = (100 * ds.counts[biotype == "rRNA"].sum() / ds.counts.sum()).round(2)
        products = ann["product"].to_dict()
        lib["top_feature"] = [
            f"{g} ({products[g]})" if g in products else g for g in ds.counts.idxmax()
        ]
        res, sf = rnaseq.run_deseq2(ds, c["factor"], c["reference"], c["test"], alpha=th["alpha"])
        own_sf = rnaseq.median_of_ratios_size_factors(ds.counts)
        galaxy = read_deseq2_table(c["galaxy_results"])
        concordance = rnaseq.compare_results(res, galaxy, th["alpha"], th["min_abs_lfc"])
        concordance["size_factor_max_abs_diff_vs_pydeseq2"] = float((own_sf - sf).abs().max())
        if c.get("galaxy_normalized"):
            gnorm = pd.read_csv(c["galaxy_normalized"], sep="\t", index_col=0)
            ours = rnaseq.normalized_counts(ds.counts, own_sf)[gnorm.columns].loc[gnorm.index]
            concordance["normalized_counts_max_rel_diff_vs_galaxy"] = float(
                ((ours - gnorm).abs() / gnorm.where(gnorm > 0)).max().max()
            )
        annotated = rnaseq.annotate(res, ann)
        sig = rnaseq.significant(annotated, th["alpha"], th["min_abs_lfc"]).sort_values("padj")
        coords, explained = rnaseq.pca(rnaseq.log_normalized(ds.counts))
        ctx.table(lib, f"rnaseq_{name}_libraries")
        ctx.table(annotated.sort_values("padj"), f"rnaseq_{name}_deseq2_all")
        ctx.table(sig, f"rnaseq_{name}_deseq2_significant")
        ctx.table(pd.Series(concordance, name="value").to_frame(), f"rnaseq_{name}_vs_galaxy")
        labels = {g: (ann.at[g, "gene"] or g) if g in ann.index else g for g in sig.index}
        figs = [
            ctx.figure(
                plotting.volcano(res, c["title"], th["alpha"], th["min_abs_lfc"], labels),
                f"rnaseq_{name}_volcano",
            ),
            ctx.figure(
                plotting.lfc_concordance(
                    res, galaxy, "PyDESeq2 vs Galaxy DESeq2", concordance["lfc_pearson"]
                ),
                f"rnaseq_{name}_concordance",
            ),
            ctx.figure(
                plotting.pca_plot(
                    coords, ds.samples[c["factor"]], explained, "Samples (top-500 variable genes)"
                ),
                f"rnaseq_{name}_pca",
            ),
        ]
        out[name] = {
            "title": c["title"],
            "libraries": lib,
            "concordance": concordance,
            "direction": rnaseq.direction_counts(res, th["alpha"], th["min_abs_lfc"]),
            "top": sig.assign(abs_lfc=sig["log2FoldChange"].abs())
            .sort_values(["padj", "abs_lfc"], ascending=[True, False])
            .head(10),
            "pca_explained": [round(float(v), 6) for v in explained],
            "figures": figs,
            "reference": c["reference"],
            "test": c["test"],
        }
    return out


def run_structure(ctx: Context) -> dict[str, Any]:
    c = ctx.cfg["structure"]
    af3 = structure.af3_summary(c["af3_dir"])
    model0 = structure.chain_ca(structure.load_structure(Path(c["af3_dir"]) / "model_0.cif"))
    query = structure.af3_query_sequence(Path(c["af3_dir"]) / "job_request.json")
    if model0.sequence != query:
        raise ValueError("AlphaFold model sequence differs from the submitted job sequence")
    pae = structure.af3_pae(Path(c["af3_dir"]) / "full_data_0.json")
    bands = structure.plddt_bands(model0.bfactors)
    ctx.table(af3, "structure_af3_models")

    a = structure.chain_ca(structure.load_structure(c["homolog_a"]))
    b = structure.chain_ca(structure.load_structure(c["homolog_b"]))
    homolog = structure.compare_chains(a, b, plddt_min=70)
    pairs, _ = structure.residue_correspondence(a, b)
    bodies = structure.rigid_domains(a, b, pairs)

    s1 = structure.chain_ca(structure.load_structure(c["spike_a"]), c["spike_chain"])
    s2 = structure.chain_ca(structure.load_structure(c["spike_b"]), c["spike_chain"])
    q1 = structure.read_seqres(c["spike_a"])[c["spike_chain"]]
    q2 = structure.read_seqres(c["spike_b"])[c["spike_chain"]]
    spike = structure.compare_chains(s1, s2, q1, q2)
    _, spike_diffs = structure.residue_correspondence(s1, s2, q1, q2)
    ctx.table(pd.Series(homolog, name="value").to_frame(), "structure_dosr_homologs")
    ctx.table(bodies, "structure_dosr_rigid_bodies", index=False)
    ctx.table(pd.Series(spike, name="value").to_frame(), "structure_spike_comparison")
    ctx.table(spike_diffs, "structure_spike_differences", index=False)
    figs = [
        ctx.figure(
            plotting.plddt_and_pae(model0.resnums, model0.bfactors, pae, c["af3_label"]),
            "structure_af3_confidence",
        ),
        ctx.figure(
            plotting.plddt_profiles(
                {
                    c["homolog_a_label"]: (a.resnums, a.bfactors),
                    c["homolog_b_label"]: (b.resnums, b.bfactors),
                },
                "AlphaFold DB per-residue confidence",
            ),
            "structure_dosr_plddt",
        ),
    ]
    rbd = spike_diffs[(spike_diffs["position"] >= 319) & (spike_diffs["position"] <= 541)]
    return {
        "af3": af3,
        "af3_bands": bands,
        "af3_len": len(model0.sequence),
        "pae_mean": round(float(pae.mean()), 2),
        "homolog": homolog,
        "bodies": bodies,
        "spike": spike,
        "spike_diffs": spike_diffs,
        "spike_rbd": rbd,
        "figures": figs,
    }


RUNNERS: dict[str, Callable[[Context], dict[str, Any]]] = {
    "qc": run_qc,
    "assembly": run_assembly,
    "annotation": run_annotation,
    "pangenome": run_pangenome,
    "rnaseq": run_rnaseq,
    "structure": run_structure,
}


def _jsonable(obj: Any) -> Any:
    if isinstance(obj, pd.DataFrame | pd.Series):
        return None
    if isinstance(obj, dict):
        return {k: v for k, v in ((k, _jsonable(v)) for k, v in obj.items()) if v is not None}
    if isinstance(obj, list | tuple):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def run_pipeline(
    cfg: dict[str, Any], stages: Iterable[str] = STAGES, write_report: bool = True
) -> dict[str, Any]:
    ctx = Context(cfg)
    results: dict[str, Any] = {}
    for stage in stages:
        if stage not in RUNNERS:
            raise ValueError(f"unknown stage {stage!r}; choose from {', '.join(STAGES)}")
        t0 = time.perf_counter()
        log.info("stage %s ...", stage)
        results[stage] = RUNNERS[stage](ctx)
        log.info("stage %s done in %.1f s", stage, time.perf_counter() - t0)
    (ctx.tables / "metrics.json").write_text(
        json.dumps(_jsonable(results), indent=2, default=str) + "\n"
    )
    if write_report and set(STAGES) <= set(results):
        from pathogenomics.report import write_report as _write

        _write(results, cfg, ctx.report_dir / "REPORT.md")
    return results
