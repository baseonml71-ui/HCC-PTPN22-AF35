"""Report archive execution limitations without executing analytical code."""
from pathlib import Path
import csv
root = Path(__file__).resolve().parents[2]
with (root / 'release/SCRIPT_ARCHIVE_INDEX.tsv').open(encoding='utf-8') as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
print(f"Archive records: {len(rows)}; NO_ARCHIVAL: {sum(r['runnable_as_is'] == 'NO_ARCHIVAL' for r in rows)}")
print('Publication provenance readiness is independent of automatic pipeline execution.')
print('Read docs/RUNNING_ARCHIVED_SCRIPTS.md before adapting any historical script.')
