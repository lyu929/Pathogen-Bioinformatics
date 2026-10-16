"""Readers for the file formats produced by Galaxy, NCBI Datasets, Roary and AlphaFold."""

from pathogenomics.io.galaxy import read_deseq2_table, read_featurecounts, read_featurecounts_dir
from pathogenomics.io.gff import GffRecord, gene_annotation_table, read_gff, sequence_regions
from pathogenomics.io.reports import (
    CutadaptReport,
    FastQCReport,
    read_cutadapt_report,
    read_fastqc_data,
)
from pathogenomics.io.seq import FastqRecord, open_text, read_fasta, read_fastq

__all__ = [
    "CutadaptReport",
    "FastQCReport",
    "FastqRecord",
    "GffRecord",
    "gene_annotation_table",
    "open_text",
    "read_cutadapt_report",
    "read_deseq2_table",
    "read_fasta",
    "read_fastq",
    "read_fastqc_data",
    "read_featurecounts",
    "read_featurecounts_dir",
    "read_gff",
    "sequence_regions",
]
