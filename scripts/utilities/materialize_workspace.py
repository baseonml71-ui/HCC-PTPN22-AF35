"""Copy only code into the ignored original module layout, without executing it."""
from pathlib import Path
import csv
import shutil

root = Path(__file__).resolve().parents[2]
workspace = root / 'data/local_workspace'
with (root / 'provenance/SCRIPT_SOURCE_BINDINGS.tsv').open(encoding='utf-8') as handle:
    bindings = list(csv.DictReader(handle, delimiter='\t'))
for row in bindings:
    source = root / row['public_file']
    target = workspace / row['source_artifact']
    assert source.resolve().is_relative_to(root)
    assert target.resolve().is_relative_to(workspace.resolve())
    if target.exists():
        assert target.read_bytes() == source.read_bytes(), f'Refusing to overwrite changed code: {target.name}'
    else:
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
print(f'Copied {len(bindings)} source files to data/local_workspace.')
print('No data downloaded and no analysis executed. Execution gaps remain unresolved.')
