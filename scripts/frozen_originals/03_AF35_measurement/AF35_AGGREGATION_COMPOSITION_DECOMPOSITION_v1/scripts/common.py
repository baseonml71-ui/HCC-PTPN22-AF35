from pathlib import Path
import datetime, hashlib, json, os
import numpy as np
import pandas as pd
ROOT=Path('.'); OUT=Path(__file__).resolve().parents[1]
A=ROOT/'AF35_COMPONENT_ABLATION_ASSOCIATION_ROBUSTNESS_v1'
B=ROOT/'AF35_MATCHED_RANDOM_CLINICAL_NULL_v1'
G=ROOT/'AF35_BROAD_IMMUNE_RESIDUAL_ASSOCIATION_v1'
H=ROOT/'AF35_RESPONSE_BLIND_PROGRAM_COHERENCE_v1'
P='AF35_AGGREGATION_COMPOSITION'; C=['ERP117672','GSE302495']
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def read(p,**kw): return pd.read_csv(p,sep='\t',float_precision='round_trip',**kw)
def save(x,n): pd.DataFrame(x).to_csv(OUT/n,sep='\t',index=False,na_rep='NA')
def write(n,s): (OUT/n).write_text(s,encoding='utf-8')
def js(n,x): write(n,json.dumps(x,ensure_ascii=False,indent=2))
def corr(x,y):
 x=x-x.mean(axis=0);y=y-y.mean(axis=0)
 return (x*y).sum(axis=0)/np.sqrt((x*x).sum(axis=0)*(y*y).sum(axis=0))
def hg(x,r):
 a=x[r];b=x[~r];n=len(x)
 return (1-3/(4*(n-2)-1))*(a.mean(axis=0)-b.mean(axis=0))/np.sqrt(((len(a)-1)*a.var(axis=0,ddof=1)+(len(b)-1)*b.var(axis=0,ddof=1))/(n-2))
def protected_paths():
 # Hash every project file except the concurrently authorized I output and dependencies.
 excluded={OUT.name,'AF35_FUNCTIONAL_ROLE_STATE_ARCHITECTURE_v1','node_modules','.pnpm','__pycache__','.git'}
 paths=[]
 for base,dirs,files in os.walk(ROOT):
  dirs[:]=[d for d in dirs if d not in excluded]
  paths.extend(Path(base)/f for f in files)
 return sorted(paths)
