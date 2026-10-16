"""Genome annotation summaries and S. aureus resistance/virulence marker screening.

Marker calls use the *protein product* text, never the gene symbol alone: RefSeq
reuses symbols for unrelated genes (e.g. ``mecA`` is also the core ClpC adaptor
protein MecA, and ``msrA`` is also methionine sulfoxide reductase), so a
symbol-only screen reports false MRSA / macrolide-resistance hits.
"""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Iterable, Mapping
from pathlib import Path

import numpy as np
import pandas as pd

from pathogenomics.io.gff import GffRecord, product_of, read_gff, sequence_regions

#: product-text patterns (case-insensitive) for clinically relevant S. aureus markers
STAPH_MARKERS: dict[str, str] = {
    "mecA (PBP2a)": r"PBP2a|PBP 2a|methicillin resistance determinant MecA",
    "mecR1": r"MecR1\b",
    "mecI": r"repressor MecI\b",
    "ccrA/ccrB": r"cassette chromosome recombinase (?:CcrA|CcrB|A|B)\b",
    "ccrC": r"cassette chromosome recombinase (?:CcrC|C)\b",
    "blaZ": r"\bBlaZ\b|beta-lactamase PC1",
    "PVL (lukS/F-PV)": r"Panton-Valentine leukocidin chain|LukF-PV|LukS-PV",
    "tst (TSST-1)": r"toxic shock syndrome toxin",
    "sea": r"enterotoxin (?:type )?A$|enterotoxin SEA\b",
    "seb": r"enterotoxin (?:type )?B$|enterotoxin SEB\b",
    "sak": r"^staphylokinase$",
    "scn": r"SCIN-A\b|^staphylococcal complement inhibitor SCIN$",
    "chp": r"chemotaxis[- ]inhibit",
    "icaA": r"\bIcaA\b|N-acetylglucosaminyltransferase IcaA",
    "cna": r"collagen adhesin",
}


def _merged_length(intervals: Iterable[tuple[str, int, int]]) -> int:
    """Total length covered by the union of ``(seqid, start, end)`` intervals (1-based, inclusive)."""
    total = 0
    by_seq: dict[str, list[tuple[int, int]]] = {}
    for seqid, start, end in intervals:
        by_seq.setdefault(seqid, []).append((start, end))
    for spans in by_seq.values():
        spans.sort()
        cur_s, cur_e = spans[0]
        for s, e in spans[1:]:
            if s <= cur_e + 1:
                cur_e = max(cur_e, e)
            else:
                total += cur_e - cur_s + 1
                cur_s, cur_e = s, e
        total += cur_e - cur_s + 1
    return total


def _is_pseudo(rec: GffRecord) -> bool:
    return (
        rec.type == "pseudogene"
        or rec.attributes.get("gene_biotype") == "pseudogene"
        or rec.attributes.get("pseudo", "").lower() == "true"
    )


def annotation_summary(
    path: str | Path, name: str | None = None, genome_length: int | None = None
) -> dict[str, object]:
    """Feature counts, coding density and annotation completeness for one GFF3 file.

    ``genome_length`` defaults to the ``##sequence-region`` total; BV-BRC exports lack
    those pragmas, so pass the assembly length explicitly for them.
    """
    records = list(read_gff(path))
    types = Counter(r.type for r in records)
    if genome_length is None:
        regions = sequence_regions(path)
        genome_length = sum(regions.values()) if regions else max(r.end for r in records)
    cds = [r for r in records if r.type == "CDS" and not _is_pseudo(r)]
    products = [product_of(r) or "" for r in cds]
    hypothetical = sum(1 for p in products if "hypothetical" in p.lower() or not p)
    gene_records = [r for r in records if r.type in {"gene", "pseudogene"}]
    pseudo = sum(1 for r in gene_records if _is_pseudo(r))
    named = sum(1 for r in cds if r.attributes.get("gene"))
    cds_lengths = np.array([r.length for r in cds]) if cds else np.array([0])
    return {
        "name": name or Path(path).name.split(".gff")[0],
        "genome_length": int(genome_length),
        "genes": len(gene_records) if gene_records else len(cds) + types["tRNA"] + types["rRNA"],
        "CDS": len(cds),
        "pseudogenes": pseudo,
        "tRNA": types["tRNA"],
        "rRNA": types["rRNA"],
        "coding_density_pct": round(
            100 * _merged_length((r.seqid, r.start, r.end) for r in cds) / genome_length, 2
        ),
        "median_cds_length_bp": int(np.median(cds_lengths)),
        "hypothetical_pct": round(100 * hypothetical / max(len(cds), 1), 1),
        "cds_with_gene_symbol_pct": round(100 * named / max(len(cds), 1), 1),
    }


def marker_matrix(
    gffs: Mapping[str, str | Path], markers: Mapping[str, str] = STAPH_MARKERS
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Presence/absence of marker products per genome.

    Returns ``(presence, evidence)``: a boolean genomes x markers table and a long table
    of the matching ``(genome, marker, locus, product)`` records for auditing.
    """
    compiled = {m: re.compile(p, flags=re.IGNORECASE) for m, p in markers.items()}
    presence = pd.DataFrame(False, index=list(gffs), columns=list(markers))
    evidence: list[dict[str, str]] = []
    for genome, path in gffs.items():
        for rec in read_gff(path, types={"CDS"}):
            product = product_of(rec)
            if not product:
                continue
            for marker, pattern in compiled.items():
                if pattern.search(product.strip()):
                    presence.loc[genome, marker] = True
                    evidence.append(
                        {
                            "genome": genome,
                            "marker": marker,
                            "locus": rec.attributes.get("locus_tag", rec.attributes.get("ID", "")),
                            "symbol": rec.attributes.get("gene", ""),
                            "product": product,
                        }
                    )
    return presence, pd.DataFrame(
        evidence, columns=["genome", "marker", "locus", "symbol", "product"]
    )


def symbol_collisions(
    gffs: Mapping[str, str | Path], symbols: Iterable[str] = ("mecA", "msrA")
) -> pd.DataFrame:
    """List CDSs whose gene *symbol* is a marker name, with their actual products.

    Documents why :func:`marker_matrix` does not trust symbols (see module docstring).
    """
    wanted = set(symbols)
    rows = []
    for genome, path in gffs.items():
        for rec in read_gff(path, types={"CDS"}):
            sym = rec.attributes.get("gene")
            if sym in wanted:
                rows.append({"genome": genome, "symbol": sym, "product": product_of(rec) or ""})
    return pd.DataFrame(rows, columns=["genome", "symbol", "product"])
