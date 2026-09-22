from pathlib import Path
import hashlib, json, datetime, shutil, os
import pandas as pd
ROOT=Path('.'); OUT=Path(__file__).resolve().parents[1]
IS_G=OUT.name.startswith('AF35_BROAD'); PREFIX='AF35_BROAD_IMMUNE' if IS_G else 'AF35_PROGRAM_COHERENCE'
A=ROOT/'AF35_COMPONENT_ABLATION_ASSOCIATION_ROBUSTNESS_v1'
B=ROOT/'AF35_MATCHED_RANDOM_CLINICAL_NULL_v1'
E=ROOT/'CLINICAL_ANCHOR_READOUT_DECOMPOSITION_v1'
P=ROOT/'HCC_ICI_project'
COHORTS=['ERP117672','GSE302495']
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def save(x,name): pd.DataFrame(x).to_csv(OUT/name,sep='\t',index=False,na_rep='NA')
def write(name,s): (OUT/name).write_text(s,encoding='utf-8')
def js(name,x): write(name,json.dumps(x,ensure_ascii=False,indent=2))
def read(p,**kw): return pd.read_csv(p,sep='\t',**kw)
def scope():
 names=[A.name,B.name,E.name,'WU_AF35_SPATIAL_CELLULAR_MICROENV_CONTEXT_v1','WU_AF35_DC_CD8_SPATIAL_TOPOLOGY_v1','AF35_FUNCTIONAL_BLOCK_ABLATION_ROBUSTNESS_v1','GSE235863_CLONE_STATE_DECOMPOSITION_v1','绘图视觉','HCC_ICI_project/06_figures','HCC_ICI_project/07_manuscript','写稿蓝图','_v20R2_source_completion','HCC_ICI_FIGURE3_TEXT_SYNC_V2_1H_2026-09-06','HCC_ICI_BIOLOGICAL_RESOLUTION_BR1_v1']
 paths=set(); excluded={'node_modules','.pnpm','__pycache__','.git'}
 for name in names:
  for base,dirs,files in os.walk(ROOT/name):
   dirs[:]=[d for d in dirs if d not in excluded]
   paths.update(Path(base)/f for f in files)
 paths.update(p for p in ROOT.iterdir() if p.is_file() and any(s in p.name.lower() for s in ['manuscript','legend','af35','checkpoint']))
 return paths,names,excluded
def initialize():
 for d in ['scripts','inputs','QA','figures','authority']: (OUT/d).mkdir(exist_ok=True)
 assert not (OUT/'QA/INITIALIZATION_LOCK.json').exists(),'Already initialized'
 pass  # Non-computational private authorization copy excluded.
 sources=[A/'AF35_ORDERED_GENES.tsv',A/'SCORING_RECOVERY_LOCK.json',B/'AF35_ORDERED_GENES.tsv',B/'MATCHED_PROGRAM_FREEZE_LOCK.json',B/'OUTCOME_JOIN_TIMESTAMP.txt',B/'COMPLETE_SHA256_MANIFEST.tsv',E/'MOLECULAR_GEOMETRY_LOCK.json']
 sources += [A/f'{c}_GENE_Z_AND_FULL_SCORES_RESPONSE_BLIND.tsv' for c in COHORTS]
 sources += [E/f'inputs/{c}_MOLECULAR_INPUT_RESPONSE_BLIND.tsv' for c in COHORTS]
 if IS_G:
  sources += [E/'SOURCE_BINDINGS.tsv',E/'scripts/02_geometry_response_blind.R',P/'01_metadata/PHASE3_PTPN22_SPECIFICITY_COMPARATOR_PROGRAMS.tsv',P/'01_metadata/phase15c_public_benchmark/PHASE15C_BENCHMARK_CONTRACT.tsv',P/'04_scripts/phase15c_public_benchmark/01_run_phase15c_benchmark.R',P/'03_processed_data/pretreatment_clinical_translation/ERP117672_gene_expression_response_blind.rds',P/'03_processed_data/phase15b_public_replication/GSE302495_RESPONSE_BLIND_MATRIX_AND_SCORES.rds']
 else:
  sources += [B/'AF35_MATCHED_NULL_PROGRAM_MANIFEST.tsv',B/'SCORING_AUTHORITY_LOCK.tsv',B/'PATIENT_SCORE_FREEZE_SHA256.tsv']+[B/f'{c}_SELECTED_TRANSFORMED_EXPRESSION.tsv' for c in COHORTS]
 assert all(p.is_file() for p in sources)
 paths,names,excluded=scope();paths.update(sources)
 save([dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(paths)],'QA/PROTECTED_BEFORE.tsv')
 js('QA/INITIALIZATION_LOCK.json',dict(utc=now(),protected_directories=names,excluded_dependency_components=sorted(excluded),protected_files=len(paths),module=PREFIX))
 # Verify direct sources against previously frozen hashes; never parse clinical tables.
 expected={}
 for base,lock,key in [(A,'SCORING_RECOVERY_LOCK.json','hashes'),(B,'MATCHED_PROGRAM_FREEZE_LOCK.json','hashes'),(E,'MOLECULAR_GEOMETRY_LOCK.json','files')]:
  expected.update({str(base/n):h for n,h in json.loads((base/lock).read_text())[key].items()})
 for base in [A,B,E]:
  m=read(base/'COMPLETE_SHA256_MANIFEST.tsv')
  expected.update({str(base/r.path):r.sha256 for r in m.itertuples()})
 if IS_G:
  expected.update({str(ROOT/r.path):r.sha256 for r in read(E/'SOURCE_BINDINGS.tsv').itertuples()})
 rows=[dict(path=str(p),sha256=sha(p),historical_sha256=expected.get(str(p),'LOCK_OR_MANIFEST_AUTHORITY'),status='PASS' if str(p) not in expected or sha(p)==expected[str(p)] else 'FAIL',content_access='HASH_ONLY_UNLESS_RESPONSE_BLIND_ALLOWLIST') for p in sources]
 save(rows,PREFIX+'_SOURCE_BINDING_AUDIT.tsv');assert all(r['status']=='PASS' for r in rows)
 g=read(A/'AF35_ORDERED_GENES.tsv');assert len(g)==35 and g.gene.nunique()==35 and 'PTPN22' not in g.gene.tolist()
 assert g.gene.tolist()==read(B/'AF35_ORDERED_GENES.tsv').gene.tolist()
 save(g,'inputs/AF35_ORDERED_GENES.tsv')
 for c in COHORTS:
  ids=read(E/f'inputs/{c}_MOLECULAR_INPUT_RESPONSE_BLIND.tsv',usecols=['patient_id']).patient_id
  assert len(ids)==(35 if c==COHORTS[0] else 38) and ids.is_unique
  save(pd.DataFrame({'patient_id':ids}),'inputs/'+c+'_PATIENT_IDS.tsv')
 print(PREFIX,'INITIALIZED',len(paths),'protected;',len(rows),'bound sources',flush=True)
if __name__=='__main__':initialize()
