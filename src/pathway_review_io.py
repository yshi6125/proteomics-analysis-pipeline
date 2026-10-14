"""Local resource validation and provenance for optional pathway review stages."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def safe_label(label: str) -> str:
    if not label or Path(label).name != label or label in {'.', '..'}:
        raise ValueError('comparison must be a filename-safe label')
    return label


def sha256_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dataframe_hash(frame: pd.DataFrame) -> str:
    return hashlib.sha256(frame.to_csv(index=False).encode()).hexdigest()


def load_definitions(path: Path) -> dict:
    raw = json.loads(Path(path).read_text())
    if not isinstance(raw, dict) or not raw:
        raise ValueError('Gene sets must be a nonempty library -> term -> gene-list JSON object')
    normalized = {}
    for library, terms in raw.items():
        if not isinstance(library, str) or not library or not isinstance(terms, dict):
            raise ValueError('Invalid library definitions')
        normalized[library] = {}
        for term, genes in terms.items():
            if not isinstance(term, str) or not term or not isinstance(genes, list):
                raise ValueError('Each term must have a gene-symbol list')
            if any(not isinstance(gene, str) or not gene.strip() for gene in genes):
                raise ValueError('Gene symbols must be nonempty strings')
            normalized[library][term] = sorted({gene.strip().upper() for gene in genes})
    return normalized


def provenance(paths) -> dict:
    return {str(Path(path).resolve()): sha256_file(path) for path in sorted(set(map(Path, paths)))}


def guard_output(output_dir: Path, source_dirs, source_files=()) -> Path:
    """Reject writes inside input directories or on top of source files."""
    target = Path(output_dir).resolve()
    for directory in source_dirs:
        source = Path(directory).resolve()
        if target == source or target.is_relative_to(source) or source.is_relative_to(target):
            raise ValueError(f'Output directory overlaps protected source directory: {source}')
    for path in source_files:
        if Path(path).resolve().is_relative_to(target):
            raise ValueError(f'Output directory contains an input file: {path}')
    return target
