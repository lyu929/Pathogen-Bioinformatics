"""Minimal, dependency-free FASTA/FASTQ readers (plain or gzip-compressed)."""

from __future__ import annotations

import gzip
import io
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

PathLike = str | Path


def open_text(path: PathLike) -> TextIO:
    """Open a text file, transparently decompressing ``.gz`` (detected by magic bytes)."""
    path = Path(path)
    with path.open("rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_fasta(path: PathLike) -> Iterator[tuple[str, str]]:
    """Yield ``(header, sequence)`` pairs. The header excludes the leading ``>``."""
    header: str | None = None
    chunks: list[str] = []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            else:
                if header is None:
                    raise ValueError(f"{path}: sequence data before first FASTA header")
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


@dataclass(frozen=True, slots=True)
class FastqRecord:
    name: str
    sequence: str
    quality: str

    def phred(self, offset: int = 33) -> list[int]:
        return [ord(c) - offset for c in self.quality]


def read_fastq(path: PathLike) -> Iterator[FastqRecord]:
    """Yield records from a 4-line FASTQ file, validating the record structure."""
    with open_text(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline().rstrip("\r\n")
            plus = fh.readline()
            qual = fh.readline().rstrip("\r\n")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"{path}: malformed FASTQ record near {header.strip()!r}")
            if len(seq) != len(qual):
                raise ValueError(f"{path}: sequence/quality length mismatch in {header.strip()!r}")
            yield FastqRecord(header[1:].strip(), seq, qual)
