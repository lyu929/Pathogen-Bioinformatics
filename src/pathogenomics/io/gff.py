"""GFF3 parsing for NCBI RefSeq, Prokka and BV-BRC (RAST) annotation files."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

from pathogenomics.io.seq import open_text


@dataclass(slots=True)
class GffRecord:
    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: dict[str, str] = field(default_factory=dict)

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def parse_attributes(column: str) -> dict[str, str]:
    """Parse a GFF3 column-9 string; values are percent-decoded (``%2C`` -> ``,``)."""
    attrs: dict[str, str] = {}
    for item in column.strip().split(";"):
        if not item:
            continue
        key, sep, value = item.partition("=")
        if not sep:
            continue
        attrs[key.strip()] = unquote(value.strip())
    return attrs


def read_gff(path: str | Path, types: Iterable[str] | None = None) -> Iterator[GffRecord]:
    """Yield feature records, optionally restricted to the given feature ``types``.

    Parsing stops at a ``##FASTA`` directive (Prokka appends sequences after it).
    """
    wanted = set(types) if types is not None else None
    with open_text(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            if line.startswith("##FASTA"):
                return
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) != 9:
                raise ValueError(
                    f"{path}:{lineno}: expected 9 tab-separated columns, got {len(cols)}"
                )
            if wanted is not None and cols[2] not in wanted:
                continue
            yield GffRecord(
                seqid=cols[0],
                source=cols[1],
                type=cols[2],
                start=int(cols[3]),
                end=int(cols[4]),
                score=cols[5],
                strand=cols[6],
                phase=cols[7],
                attributes=parse_attributes(cols[8]),
            )


def sequence_regions(path: str | Path) -> dict[str, int]:
    """Return ``{seqid: length}`` from ``##sequence-region`` pragmas (empty if absent)."""
    regions: dict[str, int] = {}
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4:
                    regions[parts[1]] = int(parts[3]) - int(parts[2]) + 1
            elif line.startswith("##FASTA"):
                break
    return regions


def product_of(record: GffRecord) -> str | None:
    """Best-effort protein product: RefSeq/Prokka ``product``; BV-BRC stores it in ``Name``."""
    if "product" in record.attributes:
        return record.attributes["product"]
    if record.source.startswith("rast") and "Name" in record.attributes:
        return record.attributes["Name"]
    return None


def gene_annotation_table(path: str | Path) -> pd.DataFrame:
    """Locus-tag keyed table of gene symbol, biotype, product and coordinates.

    Built from ``gene``/``pseudogene`` features, with the product taken from the
    child CDS/RNA feature (RefSeq stores products on the child, not the gene).
    """
    genes: dict[str, dict[str, object]] = {}
    gene_id_to_locus: dict[str, str] = {}
    products: dict[str, str] = {}
    for rec in read_gff(path):
        attrs = rec.attributes
        if rec.type in {"gene", "pseudogene"}:
            locus = attrs.get("locus_tag") or attrs.get("ID", "")
            if not locus:
                continue
            genes[locus] = {
                "locus_tag": locus,
                "gene": attrs.get("gene", ""),
                "biotype": attrs.get("gene_biotype", rec.type),
                "seqid": rec.seqid,
                "start": rec.start,
                "end": rec.end,
                "strand": rec.strand,
            }
            if "ID" in attrs:
                gene_id_to_locus[attrs["ID"]] = locus
        else:
            product = product_of(rec)
            if not product:
                continue
            locus = attrs.get("locus_tag") or gene_id_to_locus.get(attrs.get("Parent", ""), "")
            if locus and locus not in products:
                products[locus] = product
    table = pd.DataFrame.from_records(list(genes.values()))
    if table.empty:
        return pd.DataFrame(
            columns=["gene", "biotype", "product", "seqid", "start", "end", "strand"]
        )
    table["product"] = table["locus_tag"].map(products).fillna("")
    return table.set_index("locus_tag")
