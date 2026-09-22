from pathlib import Path
import hashlib,json,datetime
import pandas as pd

ROOT=Path('.')
OUT=ROOT/'WU_AF35_DC_CD8_SPATIAL_TOPOLOGY_v1'
C=ROOT/'WU_AF35_SPATIAL_CELLULAR_MICROENV_CONTEXT_v1'
DC='Dendritic cells'
CD8='CD8 T cells'
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()
def read(p):return pd.read_csv(p,sep='\t',dtype={'patient_id':str,'cell_id':str,'region_id':str,'spot_id':str})
def save(d,name):d.to_csv(OUT/name,sep='\t',index=False,na_rep='NA',float_format='%.16g')
def write_json(name,data):(OUT/name).write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def verify_lock():
    lock=json.loads((OUT/'WU_AF35_DC_CD8_TOPOLOGY_FREEZE_LOCK.json').read_text(encoding='utf-8'))
    assert sha(OUT/'WU_AF35_DC_CD8_TOPOLOGY_PROTOCOL_FREEZE_v1.md')==lock['protocol_sha256']
    for p,h in lock['bound_inputs'].items():assert sha(p)==h,p
    return lock
