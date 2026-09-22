"""Report static execution gaps; never import analytical scripts."""
from pathlib import Path
import csv

root = Path(__file__).resolve().parents[2]
with (root / 'provenance/EXECUTION_GAPS.tsv').open(encoding='utf-8') as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
for row in rows:
    print(f"{row['gap_id']}: {row['status']} — {row['description']}")
print('Public reproduction is blocked. Do not fabricate inputs or bypass frozen guards.')
raise SystemExit(2 if any(row['status'] == 'UNRESOLVED' for row in rows) else 0)
