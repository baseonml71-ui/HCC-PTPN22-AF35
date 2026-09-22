from pathlib import Path
import hashlib,json,datetime,os,shutil
import pandas as pd
ROOT=Path('.');OUT=Path(__file__).resolve().parents[1];P=ROOT/'HCC_ICI_project'
A=ROOT/'AF35_COMPONENT_ABLATION_ASSOCIATION_ROBUSTNESS_v1';H=ROOT/'AF35_RESPONSE_BLIND_PROGRAM_COHERENCE_v1';C=ROOT/'GSE235863_CLONE_STATE_DECOMPOSITION_v1'
COHORTS=['ERP117672','GSE302495']
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 with open(p,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p,**kw):return pd.read_csv(p,sep='\t',**kw)
def save(x,n):pd.DataFrame(x).to_csv(OUT/n,sep='\t',index=False,na_rep='NA')
def write(n,s):(OUT/n).write_text(s,encoding='utf-8')
def js(n,x):write(n,json.dumps(x,ensure_ascii=False,indent=2))
def jl(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def scope():
 names=[A.name,H.name,C.name,'AF35_MATCHED_RANDOM_CLINICAL_NULL_v1','CLINICAL_ANCHOR_READOUT_DECOMPOSITION_v1','WU_AF35_SPATIAL_CELLULAR_MICROENV_CONTEXT_v1','WU_AF35_DC_CD8_SPATIAL_TOPOLOGY_v1','AF35_FUNCTIONAL_BLOCK_ABLATION_ROBUSTNESS_v1','AF35_BROAD_IMMUNE_RESIDUAL_ASSOCIATION_v1','绘图视觉','HCC_ICI_BIOLOGICAL_RESOLUTION_BR1_v1','HCC_ICI_project/06_figures','HCC_ICI_project/07_manuscript','写稿蓝图','_v20R2_source_completion','HCC_ICI_FIGURE3_TEXT_SYNC_V2_1H_2026-09-06','HCC_ICI_project/05_results/tables/mechanism','HCC_ICI_project/03_processed_data/mechanism','HCC_ICI_project/05_results/phase17a_public_perturbation','HCC_ICI_project/05_results/tables/phase17a_public_perturbation','HCC_ICI_project/10_intermediate_files/phase17a_public_perturbation','HCC_ICI_project/01_metadata/phase17a_public_perturbation','HCC_ICI_project/04_scripts/phase17_dual_strengthening']
 paths=set();excluded={'node_modules','.pnpm','.git','__pycache__'}
 for n in names:
  for base,dirs,files in os.walk(ROOT/n):
   dirs[:]=[d for d in dirs if d not in excluded];paths.update(Path(base)/v for v in files)
 paths.update(p for p in ROOT.iterdir() if p.is_file() and any(v in p.name.lower() for v in ['af35','manuscript','legend','checkpoint']))
 return paths,names,excluded
if __name__=='__main__':
 for d in ['scripts','inputs','authority','QA','figures']: (OUT/d).mkdir(exist_ok=True)
 shutil.copyfile('private_dependencies/unavailable_protocol.txt',OUT/'authority/USER_REQUEST.txt')
 sources=[A/'AF35_ORDERED_GENES.tsv',H/'AF35_PROGRAM_COHERENCE_FREEZE_LOCK.json',H/'COMPLETE_SHA256_MANIFEST.tsv',C/'COMPLETE_SHA256_MANIFEST.tsv',P/'03_processed_data/phase3/GSE235863_nine_patients_scRNAseq_cd45_raw_counts.h5ad',C/'authority/03_processed_data/phase3/GSE235863_BARCODE_CLONOTYPE_STATE_PTPN22.tsv.gz',C/'authority/03_processed_data/phase19a_af35_clonotype/GSE235863_AF35_RESPONSE_BLIND_PROJECTION.npz',C/'authority/01_metadata/phase19a_af35_clonotype/PHASE19A_AF35_RESPONSE_BLIND_PROJECTION_FREEZE.json',P/'01_metadata/phase17_dual_strengthening/PHASE17_FILE_MANIFEST_SHA256.tsv']
 for c in COHORTS:sources += [H/f'inputs/{c}_GENE_Z_RESPONSE_BLIND.tsv',H/f'inputs/{c}_PATIENT_IDS.tsv',A/f'{c}_GENE_Z_AND_FULL_SCORES_RESPONSE_BLIND.tsv']
 paths,names,excluded=scope();paths.update(sources);assert all(p.is_file() for p in sources)
 if not (OUT/'QA/INITIALIZATION_LOCK.json').exists():
  save([dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(paths)],'QA/PROTECTED_BEFORE.tsv')
  js('QA/INITIALIZATION_LOCK.json',dict(utc=now(),protected_directories=names,excluded_components=sorted(excluded),protected_n=len(paths)))
 expected={}
 for base in [A,H,C]:
  m=read(base/'COMPLETE_SHA256_MANIFEST.tsv');key=next(k for k in ['path','relative_path','file'] if k in m.columns)
  expected.update({str(base/r[key]):r['sha256'] for r in m.to_dict('records')})
 expected[str(sources[4])]='a1b5ce8a7f44cc50cba2458090dc079d7548f5567d1d9ee5237b791753c8d9d4'
 rows=[dict(path=str(p),sha256=sha(p),historical_sha256=expected.get(str(p),'MANIFEST_AUTHORITY'),status='PASS' if str(p) not in expected or sha(p)==expected[str(p)] else 'FAIL') for p in sources]
 save(rows,'AF35_FUNCTIONAL_STATE_SOURCE_BINDINGS.tsv');assert all(r['status']=='PASS' for r in rows)
 shutil.copyfile(A/'AF35_ORDERED_GENES.tsv',OUT/'inputs/AF35_ORDERED_GENES.tsv')
 print('I initialized;',len(paths),'protected;',len(rows),'bound sources',flush=True)
