import json
import shutil
from pathlib import Path

import pytest

from pathogenomics.cli import main
from pathogenomics.config import ConfigError, load_config
from pathogenomics.pipeline import run_pipeline

ROOT = Path(__file__).resolve().parents[1]


def test_cli_assembly_and_trim(capsys, write_text, tmp_path):
    fa = write_text("a.fa", ">c1\nACGTACGT\n>c2\nAC\n")
    assert main(["assembly-stats", str(fa)]) == 0
    assert "n50" in capsys.readouterr().out
    fq = write_text("r.fq", "@r\nACGT\n+\nIIII\n")
    assert main(["trim", str(fq), "-m", "1", "-o", str(tmp_path / "o.fq")]) == 0
    assert json.loads(capsys.readouterr().out)["written"] == 1


def test_cli_reports_errors_instead_of_tracebacks(capsys, tmp_path):
    assert main(["assembly-stats", str(tmp_path / "missing.fa")]) == 2
    assert "error:" in capsys.readouterr().err


def test_config_lists_all_missing_inputs(tmp_path):
    cfg = (ROOT / "config/course.toml").read_text()
    (tmp_path / "config").mkdir()
    (tmp_path / "config/course.toml").write_text(cfg)
    with pytest.raises(ConfigError) as exc:
        load_config(tmp_path / "config/course.toml")
    assert "data/qc/reads/BR-Tn1_S4_L001_R1_001.fastq.gz" in str(exc.value)
    assert "data/pangenome/roary_gene_presence_absence.csv.gz" in str(exc.value)


@pytest.mark.slow
def test_full_pipeline_writes_report(tmp_path):
    cfg = load_config(ROOT / "config/course.toml")
    cfg["output"] = {"results": tmp_path / "results", "report": tmp_path / "reports"}
    results = run_pipeline(cfg)
    assert results["qc"]["trim_matches_cutadapt"]
    report = (tmp_path / "reports/REPORT.md").read_text()
    assert "reproduces the Galaxy Cutadapt 5.1 report **exactly**" in report
    assert (tmp_path / "reports/figures/dark/rnaseq_mtb_volcano.png").exists()
    metrics = json.loads((tmp_path / "results/metrics.json").read_text())
    assert metrics["pangenome"]["accessory"] == 64
    shutil.rmtree(tmp_path)
