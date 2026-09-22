from pathlib import Path
import hashlib,json,shutil
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];OLD=ROOT.parent/'PTPN22_AF35_COUPLING_ARCHITECTURE_v1'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
for d in ['01_source','02_manifest','03_lineage','04_lineage_only','05_loso','06_contributions','07_sensitivity','08_figures','09_qa','10_review','11_checkpoint','data','logs']:(ROOT/d).mkdir(exist_ok=True)
pass  # Non-computational private authorization copy excluded.
manifest=pd.read_csv(OLD/'09_qa/DELIVERABLE_MANIFEST_SHA256.tsv',sep='\t');expected=manifest.set_index('path').sha256.to_dict()
verified=[]
for r in manifest.itertuples():
    p=OLD/r.path;h=sha(p);assert h==r.sha256,r.path
    verified.append(dict(path=r.path,sha256=h,PASS=True))
pd.DataFrame(verified).to_csv(ROOT/'09_qa/PREDECESSOR_INITIAL_HASH_VERIFICATION.tsv',sep='\t',index=False)
paths=['11_checkpoint/PTPN22_AF35_COUPLING_ARCHITECTURE_CHECKPOINT.md','02_manifest/PTPN22_AF35_COUPLING_ANALYSIS_MANIFEST_V1.md','03_primary/PTPN22_AF35_PATIENT_COUPLING_DECOMPOSITION.tsv','05_state_architecture/PTPN22_AF35_STATE_CENTROIDS.tsv','05_state_architecture/PTPN22_AF35_STATE_SPECIFIC_COUPLING.tsv','07_robustness/PTPN22_AF35_COUPLING_LOPO.tsv','07_robustness/PTPN22_AF35_TECHNICAL_SENSITIVITY.tsv','09_qa/AF35_MEMBERSHIP_AUDIT.tsv','09_qa/CELL_COUNT_AUDIT.tsv','01_source/SOURCE_SHA256.tsv','01_source/RELEASED_STATE_ORDER.tsv','data/ANALYSIS_INPUT_RESPONSE_BLIND.tsv.gz','scripts/coupling_math.py','11_checkpoint/FINAL_GATES.json','09_qa/NUMERIC_QA_SUMMARY.json','09_qa/DELIVERABLE_MANIFEST_SHA256.tsv']
audit=[]
for rel in paths:
    p=OLD/rel;h=sha(p)
    row=dict(path=str(p),relative_path=rel,sha256=h,previous_sha256=expected.get(rel,'SELF_HASH_SIDECAR'),bytes=p.stat().st_size,reuse_status='HASH_VERIFIED',read_status='READ_COMPLETE')
    if '.tsv' in p.name:
        frame=pd.read_csv(p,sep='\t');row.update(rows=len(frame),columns=';'.join(frame.columns))
    else:row.update(rows=None,columns=None);p.read_text(encoding='utf-8')
    audit.append(row)
pd.DataFrame(audit).to_csv(ROOT/'01_source/SOURCE_REUSE_AUDIT.tsv',sep='\t',index=False,na_rep='NA')
pd.DataFrame(audit)[['path','bytes','sha256','previous_sha256']].to_csv(ROOT/'01_source/SOURCE_SHA256.tsv',sep='\t',index=False)
shutil.copyfile(OLD/'data/ANALYSIS_INPUT_RESPONSE_BLIND.tsv.gz',ROOT/'data/ANALYSIS_INPUT_RESPONSE_BLIND.tsv.gz')
shutil.copyfile(OLD/'scripts/coupling_math.py',ROOT/'scripts/coupling_math.py')
shutil.copyfile(OLD/'09_qa/AF35_MEMBERSHIP_AUDIT.tsv',ROOT/'09_qa/AF35_MEMBERSHIP_AUDIT.tsv')
st=pd.read_csv(OLD/'01_source/RELEASED_STATE_ORDER.tsv',sep='\t')
st['lineage']=st.released_state.str.extract(r'^(CD4|CD8)_')[0]
assert len(st)==20 and st.lineage.notna().all() and st.lineage.value_counts().to_dict()=={'CD8':11,'CD4':9}
st.to_csv(ROOT/'03_lineage/PTPN22_STATE_LINEAGE_MAP.tsv',sep='\t',index=False)
data=pd.read_csv(ROOT/'data/ANALYSIS_INPUT_RESPONSE_BLIND.tsv.gz',sep='\t')
assert len(data)==9088 and data.patient.nunique()==6 and data.barcode.is_unique and set(data.sub_cluster)==set(st.released_state)
assert set(data.timepoint)=={'pre'} and set(data.tissue)=={'T'}
assert not any(any(x in c.lower() for x in ['outcome','response','expansion','persistence']) for c in data.columns)
data['lineage']=data.sub_cluster.map(st.set_index('released_state').lineage)
counts=data.groupby(['patient','lineage']).agg(cell_n=('barcode','size'),state_n=('sub_cluster','nunique')).reset_index()
counts.to_csv(ROOT/'09_qa/CELL_COUNT_AUDIT.tsv',sep='\t',index=False)
st.assign(PASS=True).to_csv(ROOT/'09_qa/LINEAGE_MAPPING_AUDIT.tsv',sep='\t',index=False)
print('PREDECESSOR HASHES VERIFIED',len(verified));print(counts.to_string(index=False))
print('READ COMPLETED',len(paths),'sources; no new coupling estimates calculated')
