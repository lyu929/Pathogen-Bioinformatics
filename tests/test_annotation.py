import glob

import pandas as pd
import pytest

from pathogenomics.annotation import annotation_summary, marker_matrix, symbol_collisions

SYNTHETIC = """##gff-version 3
##sequence-region c 1 1000
c\tRefSeq\tgene\t1\t100\t.\t+\t.\tID=g1;gene=mecA;locus_tag=L1;gene_biotype=protein_coding
c\tRefSeq\tCDS\t1\t100\t.\t+\t0\tParent=g1;gene=mecA;locus_tag=L1;product=adaptor protein MecA
c\tRefSeq\tgene\t51\t200\t.\t+\t.\tID=g2;locus_tag=L2;gene_biotype=protein_coding
c\tRefSeq\tCDS\t51\t200\t.\t+\t0\tParent=g2;locus_tag=L2;product=hypothetical protein
c\tRefSeq\tgene\t301\t400\t.\t+\t.\tID=g3;locus_tag=L3;gene_biotype=pseudogene
c\tRefSeq\tCDS\t301\t400\t.\t+\t0\tParent=g3;locus_tag=L3;pseudo=true;product=staphylokinase
c\tRefSeq\ttRNA\t500\t570\t.\t+\t.\tlocus_tag=L4;product=tRNA-Gly
"""


def test_annotation_summary_counts_and_density(write_text):
    s = annotation_summary(write_text("s.gff", SYNTHETIC), "toy")
    assert s["genome_length"] == 1000
    assert (s["genes"], s["CDS"], s["pseudogenes"], s["tRNA"]) == (3, 2, 1, 1)
    # union of CDS 1-100 and 51-200 = 200 bp of 1000; the pseudogene CDS is excluded
    assert s["coding_density_pct"] == pytest.approx(20.0)
    assert s["hypothetical_pct"] == pytest.approx(50.0)


def test_symbol_alone_never_calls_a_marker(write_text):
    path = write_text("s.gff", SYNTHETIC)
    presence, evidence = marker_matrix({"toy": path})
    assert not presence.loc["toy", "mecA (PBP2a)"]
    assert presence.loc["toy", "sak"]  # matched on product text (pseudo CDS still annotated)
    assert symbol_collisions({"toy": path}).to_dict("records") == [
        {"genome": "toy", "symbol": "mecA", "product": "adaptor protein MecA"}
    ]
    assert set(evidence["marker"]) == {"sak"}


def test_bvbrc_products_live_in_name(write_text):
    gff = write_text(
        "b.gff",
        "Unknown\trast2_export\tCDS\t1\t90\t.\t+\t0\tID=fig|1.peg.1;Name=Penicillin-binding protein PBP2a, methicillin resistance determinant MecA\n",
    )
    presence, _ = marker_matrix({"b": gff})
    assert presence.loc["b", "mecA (PBP2a)"]
    assert annotation_summary(gff, genome_length=100)["coding_density_pct"] == pytest.approx(90.0)


def test_course_staph_mrsa_calls(data_dir):
    strains = pd.read_csv(
        data_dir / "annotation/staphylococcus_aureus/refseq/strains.tsv", sep="\t"
    ).set_index("accession")
    gffs = {
        strains.loc[p.split("/")[-1].split(".gff")[0], "strain"]: p
        for p in sorted(
            glob.glob(str(data_dir / "annotation/staphylococcus_aureus/refseq/*.gff.gz"))
        )
    }
    presence, _ = marker_matrix(gffs)
    mrsa = set(presence.index[presence["mecA (PBP2a)"]])
    assert mrsa == {"MRSA107", "M48", "0073_III_ST239", "HU-14"}
    # FDAARGOS_773 carries a gene *named* mecA, but it is the ClpC adaptor, not PBP2a
    collisions = symbol_collisions({"FDAARGOS_773": gffs["FDAARGOS_773"]}, symbols=["mecA"])
    assert collisions["product"].tolist() == ["adaptor protein MecA"]
