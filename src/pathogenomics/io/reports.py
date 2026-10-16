"""Parsers for FastQC ``fastqc_data.txt`` and Cutadapt text reports."""

from __future__ import annotations

import io
import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from pathogenomics.io.seq import open_text


@dataclass
class FastQCReport:
    """Parsed ``fastqc_data.txt``: module statuses, basic statistics and module tables."""

    version: str
    statuses: dict[str, str]
    basic: dict[str, str]
    tables: dict[str, pd.DataFrame] = field(default_factory=dict)

    @property
    def total_sequences(self) -> int:
        return int(self.basic["Total Sequences"])

    @property
    def gc_percent(self) -> int:
        return int(self.basic["%GC"])

    @property
    def sequence_length(self) -> str:
        return self.basic["Sequence length"]

    def per_base_quality(self) -> pd.DataFrame:
        """Per-base quality table with ``start``/``end`` columns parsed from FastQC bins."""
        table = self.tables["Per base sequence quality"].copy()
        bounds = table["Base"].astype(str).str.split("-", expand=True)
        table["start"] = bounds[0].astype(int)
        table["end"] = bounds[bounds.columns[-1]].fillna(bounds[0]).astype(int)
        return table


def read_fastqc_data(path: str | Path) -> FastQCReport:
    version = ""
    statuses: dict[str, str] = {}
    tables: dict[str, pd.DataFrame] = {}
    module: str | None = None
    buffer: list[str] = []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith("##FastQC"):
                version = line.split("\t")[-1]
            elif line.startswith(">>END_MODULE"):
                if module is not None:
                    header = next((b for b in buffer if b.startswith("#")), None)
                    rows = [b for b in buffer if not b.startswith("#")]
                    if header is not None and rows:
                        text = "\n".join([header.lstrip("#"), *rows])
                        tables[module] = pd.read_csv(io.StringIO(text), sep="\t")
                module, buffer = None, []
            elif line.startswith(">>"):
                name, _, status = line[2:].partition("\t")
                module, buffer = name, []
                statuses[name] = status
            elif module is not None:
                # keep only the column header ("#Base ...") and data rows; drop extra comments
                if line.startswith("#") and any(b.startswith("#") for b in buffer):
                    continue
                buffer.append(line)
    if "Basic Statistics" not in tables:
        raise ValueError(f"{path}: no 'Basic Statistics' module - not a fastqc_data.txt file")
    basic_table = tables.pop("Basic Statistics")
    basic = dict(
        zip(basic_table.iloc[:, 0].astype(str), basic_table.iloc[:, 1].astype(str), strict=True)
    )
    return FastQCReport(version=version, statuses=statuses, basic=basic, tables=tables)


@dataclass(frozen=True)
class CutadaptReport:
    version: str
    command_line: str
    paired: bool
    reads_processed: int
    reads_with_adapters: int | None
    too_short: int
    written: int
    bp_processed: int
    bp_quality_trimmed: int
    bp_written: int

    @property
    def fraction_reads_kept(self) -> float:
        return self.written / self.reads_processed

    @property
    def fraction_bp_kept(self) -> float:
        return self.bp_written / self.bp_processed

    def option(self, name: str) -> str | None:
        """Value of a ``--name=value`` / ``-x value`` option on the recorded command line."""
        match = re.search(rf"(?:^|\s){re.escape(name)}(?:=|\s+)(\S+)", self.command_line)
        return match.group(1) if match else None


_INT = r"([\d,]+)"


def _grab(pattern: str, text: str, required: bool = True) -> int | None:
    match = re.search(pattern, text, flags=re.MULTILINE)
    if match is None:
        if required:
            raise ValueError(f"cutadapt report: pattern not found: {pattern}")
        return None
    return int(match.group(1).replace(",", ""))


def read_cutadapt_report(path: str | Path) -> CutadaptReport:
    text = Path(path).read_text(encoding="utf-8")
    version_match = re.search(r"This is cutadapt (\S+)", text)
    cmd_match = re.search(r"Command line parameters: (.*)", text)
    if version_match is None or cmd_match is None:
        raise ValueError(f"{path}: not a cutadapt report")
    paired = "Total read pairs processed" in text
    unit = "Pairs" if paired else "Reads"
    return CutadaptReport(
        version=version_match.group(1),
        command_line=cmd_match.group(1).strip(),
        paired=paired,
        reads_processed=_grab(rf"Total (?:read pairs|reads) processed:\s+{_INT}", text),
        reads_with_adapters=_grab(rf"Reads with adapters:\s+{_INT}", text, required=False),
        too_short=_grab(rf"{unit} that were too short:\s+{_INT}", text),
        written=_grab(rf"{unit} written \(passing filters\):\s+{_INT}", text),
        bp_processed=_grab(rf"Total basepairs processed:\s+{_INT} bp", text),
        bp_quality_trimmed=_grab(rf"Quality-trimmed:\s+{_INT} bp", text),
        bp_written=_grab(rf"Total written \(filtered\):\s+{_INT} bp", text),
    )
