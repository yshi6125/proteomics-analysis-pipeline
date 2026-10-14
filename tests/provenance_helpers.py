"""Verify immutable validation records after an explicitly authorized canonical promotion."""
import hashlib
import json
from pathlib import Path
from src.pathway_review_io import ROOT

PREVIOUS = ROOT/'archive/route_b_10k_canonical_20261014'


def assert_recorded_sources(sources):
    records=json.loads((PREVIOUS/'MANIFEST.json').read_text())['files']
    for original,expected in sources.items():
        current=Path(original)
        if current.exists() and hashlib.sha256(current.read_bytes()).hexdigest()==expected:
            continue
        matches=[r for r in records if (ROOT/r['original_path']).resolve()==current.resolve() and r['sha256']==expected]
        assert matches, f'No preserved source matching {original}'
        assert hashlib.sha256((ROOT/matches[0]['archived_path']).read_bytes()).hexdigest()==expected
