import pytest

from pathogenomics.assembly import assembly_stats, assembly_table, contig_coverage, nx


def test_nx_on_textbook_example():
    lengths = [2, 10, 4, 8, 6]  # total 30
    assert nx(lengths, 0.5) == (8, 2)  # 10 + 8 = 18 >= 15
    assert nx(lengths, 0.9) == (4, 4)  # 10 + 8 + 6 + 4 = 28 >= 27
    assert nx(lengths, 1.0) == (2, 5)


def test_ng50_uses_genome_size_and_can_be_unreached():
    assert nx([10, 8, 6], 0.5, genome_size=40) == (6, 3)  # 10 + 8 < 20 <= 24
    assert nx([10, 8, 6], 0.5, genome_size=100) == (0, 0)
    assert nx([], 0.5) == (0, 0)
    with pytest.raises(ValueError):
        nx([1], 0)


def test_contig_coverage_from_headers():
    assert contig_coverage("NODE_1_length_936619_cov_24.799470") == pytest.approx(24.79947)
    assert contig_coverage("contig00001 len=536032 cov=58.2 corr=0") == pytest.approx(58.2)
    assert contig_coverage("chr1") is None


def test_assembly_stats_small(write_text):
    fa = write_text(
        "a.fasta",
        ">NODE_1_length_8_cov_10.0\nGGGGCCCC\n>NODE_2_length_4_cov_2.0\nAATT\n>NODE_3_length_2_cov_1\nNN\n",
    )
    s = assembly_stats(fa, min_length=3)
    assert (s.contigs, s.total_length, s.largest_contig, s.n50, s.l50) == (2, 12, 8, 8, 1)
    assert s.gc_percent == pytest.approx(66.67)
    assert s.mean_coverage == pytest.approx((10 * 8 + 2 * 4) / 12, abs=0.05)
    with pytest.raises(ValueError, match="no contigs"):
        assembly_stats(fa, min_length=100)


def test_course_pseudomonas_assemblies(data_dir):
    table = assembly_table(
        sorted((data_dir / "assembly/pseudomonas_aeruginosa").glob("*.fasta.gz")), min_length=500
    )
    first = table.loc["SRR12131093"]
    assert first["n50"] == 363_515 and first["largest_contig"] == 936_619 and first["l50"] == 6
    assert table["total_length"].between(6.49e6, 6.51e6).all()
    assert table["gc_percent"].between(66.0, 66.6).all()  # P. aeruginosa ~66% GC
