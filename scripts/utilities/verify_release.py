"""Verify the release payload without importing or running analytical modules."""
from pathlib import Path
import csv
import hashlib
import json

root = Path(__file__).resolve().parents[2]
digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
checks = {}
for line in (root / 'release/SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
    expected, name = line.split('  ', 1)
    path = root / name
    assert path.resolve().is_relative_to(root), name
    assert path.is_file() and digest(path) == expected, name
    checks[name] = expected
payload = {p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file()
           and '.git' not in p.relative_to(root).parts
           and not p.relative_to(root).as_posix().startswith('data/local_workspace/')
           and '__pycache__' not in p.parts}
assert payload == set(checks) | {'release/SHA256SUMS.txt'}, 'Unlisted or missing payload'
with (root / 'release/RELEASE_MANIFEST.tsv').open(encoding='utf-8', newline='') as handle:
    manifest = list(csv.DictReader(handle, delimiter='\t'))
assert {r['relative_path'] for r in manifest} == set(checks) - {'release/RELEASE_MANIFEST.tsv'}
for row in manifest:
    path = root / row['relative_path']
    assert path.stat().st_size == int(row['size_bytes'])
    assert checks[row['relative_path']] == row['sha256']
genes = (root / 'config/AF35_genes.txt').read_text(encoding='utf-8').splitlines()
assert len(genes) == len(set(genes)) == 35 and 'PTPN22' not in genes
citation = json.loads((root / 'CITATION.cff').read_text(encoding='utf-8'))
assert citation['version'] == '1.0.0' and len(citation['authors']) == 6
with (root / 'provenance/CURRENT_ANALYSIS_ENTRYPOINTS.tsv').open(encoding='utf-8') as handle:
    entries = list(csv.DictReader(handle, delimiter='\t'))
for row in entries:
    assert row['status'] == 'ACTIVE_CURRENT'
    assert (root / row['active_script']).is_file()
    assert (root / row['environment']).is_file()
print(f'PASS: {len(payload)} payload files; {len(entries)} source entries; AF35 35/35')
print('Integrity verification only; no biological analysis or public release performed.')
