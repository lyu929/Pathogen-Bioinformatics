import pandas as pd
import pytest

from pathogenomics.io.galaxy import read_deseq2_table, read_featurecounts, read_featurecounts_dir
from pathogenomics.io.gff import gene_annotation_table, parse_attributes, read_gff, sequence_regions
from pathogenomics.io.reports import read_cutadapt_report, read_fastqc_data
from pathogenomics.io.seq import read_fasta, read_fastq

GFF = """##gff-version 3
##sequence-region chr1 1 1000
chr1\tRefSeq\tgene\t10\t300\t.\t+\t.\tID=gene-A;gene=mecA;gene_biotype=protein_coding;locus_tag=A_0001
chr1\tProtein Homology\tCDS\t10\t300\t.\t+\t0\tID=cds-1;Parent=gene-A;product=PBP2a family transpeptidase%2C MecA;locus_tag=A_0001
chr1\tRefSeq\tgene\t400\t500\t.\t-\t.\tID=gene-B;gene_biotype=tRNA;locus_tag=A_0002
chr1\ttRNAscan-SE\ttRNA\t400\t500\t.\t-\t.\tID=rna-1;Parent=gene-B;product=tRNA-Arg
##FASTA
>chr1
ACGT
"""


def test_read_fasta_multiline_and_gzip(write_text):
    path = write_text("x.fa.gz", ">c1 desc\nACG\nTT\n\n>c2\nGG\n")
    assert list(read_fasta(path)) == [("c1 desc", "ACGTT"), ("c2", "GG")]


def test_read_fasta_rejects_sequence_before_header(write_text):
    with pytest.raises(ValueError, match="before first FASTA header"):
        list(read_fasta(write_text("bad.fa", "ACGT\n>c1\nA\n")))


def test_read_fastq_validates_records(write_text):
    ok = write_text("r.fq", "@r1\nACGT\n+\nIIII\n")
    rec = next(read_fastq(ok))
    assert rec.name == "r1" and rec.phred() == [40, 40, 40, 40]
    with pytest.raises(ValueError, match="length mismatch"):
        list(read_fastq(write_text("bad.fq", "@r1\nACGT\n+\nIII\n")))
    with pytest.raises(ValueError, match="malformed"):
        list(read_fastq(write_text("bad2.fq", "r1\nACGT\n+\nIIII\n")))


def test_parse_attributes_percent_decodes():
    assert parse_attributes("ID=x;product=a%2C b;note=50%25;flag") == {
        "ID": "x",
        "product": "a, b",
        "note": "50%",
    }


def test_read_gff_stops_at_fasta_and_filters(write_text):
    path = write_text("a.gff", GFF)
    assert [r.type for r in read_gff(path)] == ["gene", "CDS", "gene", "tRNA"]
    assert [r.length for r in read_gff(path, types={"CDS"})] == [291]
    assert sequence_regions(path) == {"chr1": 1000}


def test_read_gff_rejects_wrong_column_count(write_text):
    with pytest.raises(ValueError, match="9 tab-separated"):
        list(read_gff(write_text("b.gff", "chr1\tsrc\tgene\t1\t2\n")))


def test_gene_annotation_table_takes_product_from_child(write_text):
    table = gene_annotation_table(write_text("a.gff", GFF))
    assert table.loc["A_0001", "gene"] == "mecA"
    assert table.loc["A_0001", "product"] == "PBP2a family transpeptidase, MecA"
    assert table.loc["A_0002", "product"] == "tRNA-Arg"
    assert table.loc["A_0002", "biotype"] == "tRNA"


def test_read_featurecounts_infers_sample_and_validates(write_text):
    path = write_text(
        "SRX123-[featureCounts].tsv", "Geneid\t71 Bowtie2 on dataset 85: alignments\ng1\t5\ng2\t0\n"
    )
    s = read_featurecounts(path)
    assert s.name == "SRX123" and s.to_dict() == {"g1": 5, "g2": 0}
    with pytest.raises(ValueError, match="duplicated"):
        read_featurecounts(write_text("dup.tsv", "Geneid\tc\ng1\t1\ng1\t2\n"))
    with pytest.raises(ValueError, match="negative"):
        read_featurecounts(write_text("neg.tsv", "Geneid\tc\ng1\t-1\n"))


def test_read_featurecounts_dir_requires_same_genes(tmp_path):
    (tmp_path / "SRX1.tsv").write_text("Geneid\tc\ng1\t1\ng2\t2\n")
    (tmp_path / "SRX2.tsv").write_text("Geneid\tc\ng1\t3\ng2\t4\n")
    m = read_featurecounts_dir(tmp_path)
    assert list(m.columns) == ["SRX1", "SRX2"] and m.loc["g2", "SRX2"] == 4
    (tmp_path / "SRX3.tsv").write_text("Geneid\tc\ng1\t3\ng9\t4\n")
    with pytest.raises(ValueError, match="differs"):
        read_featurecounts_dir(tmp_path)


def test_read_deseq2_table_with_and_without_header(write_text):
    row = "g1\t10.5\t1.2\t0.3\t4\t0.001\t0.01\n"
    a = read_deseq2_table(write_text("a.tsv", row))
    b = read_deseq2_table(
        write_text("b.tsv", "gene_id\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n" + row)
    )
    pd.testing.assert_frame_equal(a, b)
    assert a.loc["g1", "padj"] == pytest.approx(0.01)


def test_course_fastqc_and_cutadapt_reports(data_dir):
    rep = read_fastqc_data(data_dir / "qc/fastqc/BR-Tn1_R1.fastqc_data.txt.gz")
    assert rep.total_sequences == 5190 and rep.gc_percent == 46 and rep.version == "0.12.1"
    assert rep.statuses["Per base sequence quality"] == "fail"
    bins = rep.per_base_quality()
    assert bins.loc[bins["Base"] == "10-11", ["start", "end"]].values.tolist() == [[10, 11]]

    paired = read_cutadapt_report(data_dir / "qc/cutadapt/SRR29263899.cutadapt.txt")
    assert paired.paired and paired.reads_processed == 1_437_637 and paired.written == 1_436_805
    assert paired.option("--quality-cutoff") == "20"
    adapter = read_cutadapt_report(data_dir / "qc/cutadapt/Rough_SCFM_3.cutadapt.txt")
    assert adapter.reads_with_adapters == 765_043
    assert adapter.option("-a") == "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
