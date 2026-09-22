from pathlib import Path
import datetime, hashlib, json, os, sys, re
import numpy as np
import pandas as pd
from scipy.stats import rankdata
ROOT=Path('.');OUT=Path(__file__).resolve().parents[1]
I=ROOT/'AF35_FUNCTIONAL_ROLE_STATE_ARCHITECTURE_v1';H=ROOT/'HCC_ICI_project'
V=ROOT/'绘图视觉/HCC_ICI_FIGURE4_V2_2_LOW_INTERRUPTION_RECONSTRUCTION_v1'
P='PTPN22_AF35';Q=OUT/'QA';A=OUT/'inputs'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 with open(p,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p,**kw):return pd.read_csv(p,sep='\t',float_precision='round_trip',**kw)
def save(x,p):pd.DataFrame(x).to_csv(p,sep='\t',index=False,na_rep='NA')
def write(p,s):Path(p).write_text(s,encoding='utf-8')
def js(p,x):write(p,json.dumps(x,ensure_ascii=False,indent=2))
def jl(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def paths():
 r=[]
 for base,dirs,files in os.walk(ROOT):
  dirs[:]=[d for d in dirs if d not in {OUT.name,'__pycache__','node_modules','.pnpm','.git'}]
  r.extend(Path(base)/f for f in files)
 return sorted(r)
