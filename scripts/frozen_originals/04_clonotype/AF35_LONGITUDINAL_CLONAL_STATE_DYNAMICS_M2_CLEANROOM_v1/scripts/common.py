from pathlib import Path
import datetime, hashlib, json, os, re, sys
import numpy as np
import pandas as pd
from pandas.io.formats import csvs
from scipy.stats import rankdata
ROOT=Path('.'); OUT=Path(__file__).resolve().parents[1]
A=OUT/'01_phaseA_outcome_free_inputs'; F=OUT/'02_phaseA_frozen_outputs'
B=OUT/'03_phaseB_outcome_inputs'; C=OUT/'04_phaseB_results'; Q=OUT/'QA'
K=ROOT/'AF35_CLINICAL_STATE_OCCUPANCY_DECOMPOSITION_v1'
L=ROOT/'AF35_CLONOTYPE_NEUTRALIZED_STATE_OCCUPANCY_v1'
I=ROOT/'AF35_FUNCTIONAL_ROLE_STATE_ARCHITECTURE_v1'
D=ROOT/'GSE235863_CLONE_STATE_DECOMPOSITION_v1'
P='AF35_M2'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 with open(p,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p,**kw):return pd.read_csv(p,sep='\t',float_precision='round_trip',**kw)
def save(x,p):pd.DataFrame(x).to_csv(p,sep='\t',index=False,na_rep='NA')
def write(p,s):Path(p).write_text(s,encoding='utf-8')
def js(p,x):write(p,json.dumps(x,ensure_ascii=False,indent=2))
def jl(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def project_files():
 pp=[]
 for base,dirs,files in os.walk(ROOT):
  dirs[:]=[d for d in dirs if d not in {OUT.name,'__pycache__','node_modules','.pnpm','.git'}]
  pp.extend(Path(base)/f for f in files)
 return sorted(pp)
def sanitize_table(d,name):
 bad=[str(c) for c in d.columns if re.search(r'response|responder|nonresponder|recist|benefit|outcome|clinical_group|r_nr|label_response',str(c),re.I) or re.search(r'(^|[^a-z])(pfs|os)($|[^a-z])',str(c),re.I)]
 values=[]
 for c in d.select_dtypes(include=['object','string']).columns:
  hit=d[c].astype(str).str.strip().str.lower().isin(['r','nr','responder','nonresponder'])
  if hit.any():values.append(str(c))
 return dict(file=name,rows=len(d),columns=len(d.columns),forbidden_columns=';'.join(bad),forbidden_value_columns=';'.join(values),status='FAIL' if bad or values else 'PASS')
def install_guard():
 # Libraries are loaded before the guard. All subsequent molecular data access
 # is limited to the two Phase-A directories. Deny subprocess and network use.
 events=[]
 def hook(event,args):
  if event=='open' and isinstance(args[0],(str,bytes,os.PathLike)):
   path=Path(os.fsdecode(args[0])).resolve(); mode=str(args[1]); flags=args[2] or 0
   writing=any(c in mode for c in 'wax+') or bool(flags & (os.O_WRONLY|os.O_RDWR|os.O_CREAT))
   allowed=path.is_relative_to(F) if writing else (path.is_relative_to(A) or path.is_relative_to(F))
   events.append(dict(path=str(path),writing=writing,allowed=allowed))
   if not allowed:raise PermissionError('M2_PHASEA_READ_WRITE_ALLOWLIST_DENIED: '+str(path))
  elif event.startswith(('subprocess.','socket.')) or event in ['os.system','os.spawn']:
   raise PermissionError('M2_PHASEA_PROCESS_NETWORK_DENIED')
 sys.addaudithook(hook)
 return events
def verify_phasea():
 lock=jl(OUT/(P+'_PHASEA_FREEZE_LOCK.json'))
 for n,h in lock['hashes'].items():assert sha(OUT/n)==h, 'FAIL_POST_FREEZE_INTEGRITY '+n
 return lock
