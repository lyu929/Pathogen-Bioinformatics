"""pathogenomics: downstream analysis of bacterial-pathogen sequencing data.

The package re-implements, in tested Python, the analyses performed in Galaxy
for GSU BIOL 4640 (Fundamentals of Bioinformatics, Fall 2025):

* read QC and Cutadapt-compatible quality trimming        -> :mod:`pathogenomics.qc`
* de novo assembly statistics                              -> :mod:`pathogenomics.assembly`
* genome annotation summaries and marker-gene screening    -> :mod:`pathogenomics.annotation`
* Roary pangenome partitioning and rarefaction             -> :mod:`pathogenomics.pangenome`
* RNA-seq differential expression (PyDESeq2) + validation  -> :mod:`pathogenomics.rnaseq`
* AlphaFold / PDB structure confidence and superposition   -> :mod:`pathogenomics.structure`
"""

__version__ = "1.0.0"
