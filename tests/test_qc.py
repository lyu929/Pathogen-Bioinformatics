import numpy as np
import pytest

from pathogenomics.io.reports import read_cutadapt_report, read_fastqc_data
from pathogenomics.qc import (
    binned_mean_quality,
    cutadapt_summary,
    fastq_stats,
    fastqc_summary,
    flag_failing_modules,
    per_base_quality_profile,
    quality_trim_index,
    trim_fastq,
)


def q(*scores: int) -> str:
    return "".join(chr(s + 33) for s in scores)


def test_trim_index_matches_cutadapt_documentation_example():
    # Cutadapt user guide, "Quality trimming algorithm": cutoff 10 keeps the first four bases
    assert quality_trim_index(q(42, 40, 26, 27, 8, 7, 11, 4, 2, 3), 10) == (0, 4)


def test_trim_index_keeps_good_reads_and_drops_hopeless_ones():
    assert quality_trim_index(q(30, 30, 30), 20) == (0, 3)
    assert quality_trim_index(q(2, 2, 2), 20) == (0, 0)
    assert quality_trim_index("", 20) == (0, 0)


def test_trim_index_front_end():
    assert quality_trim_index(q(2, 3, 30, 30, 30), 0, cutoff_front=10) == (2, 5)
    assert quality_trim_index(q(2, 30, 30, 2), 10, cutoff_front=10) == (1, 3)


def test_trim_fastq_counts_and_output(write_text, tmp_path):
    fq = write_text(
        "r.fq",
        f"@a\nACGTACGT\n+\n{q(30, 30, 30, 30, 30, 30, 2, 2)}\n@b\nACGT\n+\n{q(2, 2, 2, 2)}\n",
    )
    out = tmp_path / "o.fq"
    s = trim_fastq(fq, quality_cutoff=20, min_length=3, output=out)
    assert (s.reads_processed, s.too_short, s.written) == (2, 1, 1)
    assert (s.bp_processed, s.bp_quality_trimmed, s.bp_written) == (12, 6, 6)
    assert out.read_text() == f"@a\nACGTAC\n+\n{q(30, 30, 30, 30, 30, 30)}\n"


def test_trimming_reproduces_galaxy_cutadapt_report_exactly(data_dir):
    """Re-trimming the course FASTQ gives Cutadapt 5.1's numbers base for base."""
    report = read_cutadapt_report(data_dir / "qc/cutadapt/BR-Tn1.cutadapt.txt")
    summary = trim_fastq(
        data_dir / "qc/reads/BR-Tn1_S4_L001_R1_001.fastq.gz", quality_cutoff=20, min_length=22
    )
    assert summary.matches(report)
    assert summary.bp_quality_trimmed == 17_416


def test_fastq_stats_agree_with_fastqc(data_dir):
    stats = fastq_stats(data_dir / "qc/reads/BR-Tn1_S4_L001_R1_001.fastq.gz")
    fastqc = read_fastqc_data(data_dir / "qc/fastqc/BR-Tn1_R1.fastqc_data.txt.gz")
    assert stats.reads == fastqc.total_sequences
    assert int(stats.gc_percent) == fastqc.gc_percent  # FastQC truncates
    bins = fastqc.per_base_quality()
    np.testing.assert_allclose(binned_mean_quality(stats, bins), bins["Mean"].to_numpy(), atol=1e-9)


def test_per_base_profile_handles_ragged_reads():
    np.testing.assert_allclose(per_base_quality_profile([q(10, 20), q(30)]), [20.0, 20.0])


def test_summary_tables(data_dir):
    reports = {
        p.name.split(".")[0]: read_fastqc_data(p)
        for p in sorted((data_dir / "qc/fastqc").glob("*.gz"))
    }
    table = fastqc_summary(reports)
    fails = flag_failing_modules(table)
    assert "Adapter Content" in fails["BR-Tn1_R1"]
    assert fails["SRR10828193_R1"] == ""
    cut = cutadapt_summary(
        {p.stem: read_cutadapt_report(p) for p in (data_dir / "qc/cutadapt").glob("*.txt")}
    )
    assert cut.loc["Rough_SCFM_3.cutadapt", "pct_with_adapter"] == pytest.approx(76.5)
