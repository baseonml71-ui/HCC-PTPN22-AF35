from pathlib import Path
import datetime, hashlib, json, os
import numpy as np
import pandas as pd
ROOT=Path('.');OUT=Path(__file__).resolve().parents[1]
I=ROOT/'AF35_FUNCTIONAL_ROLE_STATE_ARCHITECTURE_v1';D=ROOT/'GSE235863_CLONE_STATE_DECOMPOSITION_v1'
P='AF35_STATE_OCCUPANCY'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 with open(p,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p,**kw):return pd.read_csv(p,sep='\t',float_precision='round_trip',**kw)
def save(x,n):pd.DataFrame(x).to_csv(OUT/n,sep='\t',index=False,na_rep='NA')
def write(n,s):(OUT/n).write_text(s,encoding='utf-8')
def js(n,x):write(n,json.dumps(x,ensure_ascii=False,indent=2))
def jl(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def paths():
 out=[]
 for base,dirs,files in os.walk(ROOT):
  dirs[:]=[d for d in dirs if d not in {OUT.name,'__pycache__','node_modules','.pnpm','.git'}]
  out.extend(Path(base)/n for n in files)
 return sorted(out)
def model(pb,patients,states):
 observed_states=[s for s in states if s in set(pb.released_state)]
 n=len(patients);s=len(observed_states);si={v:i for i,v in enumerate(observed_states)};pi={v:i for i,v in enumerate(patients)}
 # Sum-to-zero contrast coding; each observed patient-state is one unweighted row.
 contrast_s=np.vstack([np.eye(s-1),-np.ones((1,s-1))])
 contrast_p=np.vstack([np.eye(n-1),-np.ones((1,n-1))])
 mat=np.column_stack([np.ones(len(pb)),contrast_s[[si[v] for v in pb.released_state]],contrast_p[[pi[v] for v in pb.patient]]])
 adjacent={('p',v):set() for v in patients};adjacent.update({('s',v):set() for v in observed_states})
 for r in pb.itertuples():adjacent[('p',r.patient)].add(('s',r.released_state));adjacent[('s',r.released_state)].add(('p',r.patient))
 unseen=set(adjacent);components=[]
 while unseen:
  stack=[min(unseen)];vis=set()
  while stack:
   v=stack.pop()
   if v in vis:continue
   vis.add(v);stack.extend(adjacent[v]-vis)
  unseen-=vis;components.append(sorted(vis))
 rank=int(np.linalg.matrix_rank(mat));cond=float(np.linalg.cond(mat));ident=len(components)==1 and rank==mat.shape[1]
 info=dict(patient_n=n,state_n=s,observed_patient_state_n=len(pb),design_columns=mat.shape[1],rank=rank,condition_number=cond,connected_components=len(components),component_members=json.dumps(components),identifiable=ident)
 if not ident:return info,None
 coef=np.linalg.lstsq(mat,pb.m_is.to_numpy(),rcond=None)[0]
 alpha=contrast_s@coef[1:s];gamma=contrast_p@coef[s:];eps=pb.m_is.to_numpy()-mat@coef
 residuals=pb.copy();residuals['mu']=coef[0];residuals['alpha_s']=residuals.released_state.map(dict(zip(observed_states,alpha)));residuals['gamma_i']=residuals.patient.map(dict(zip(patients,gamma)));residuals['epsilon_is']=eps;residuals['fitted_m_is']=mat@coef
 rows=[]
 for pat in patients:
  d=residuals[residuals.patient==pat];m=float((d.p_is*d.m_is).sum());o=float((d.p_is*d.alpha_s).sum());w=float(d.gamma_i.iloc[0]);x=float((d.p_is*d.epsilon_is).sum())
  rows.append(dict(patient=pat,compartment=pb.compartment.iloc[0],cell_n=int(d.cell_n.sum()),observed_state_n=len(d),mu=coef[0],M_i=m,O_i=o,W_i=w,X_i=x,reconstruction_error=m-(coef[0]+o+w+x)))
 return info,dict(patient_components=pd.DataFrame(rows),state_effects=pd.DataFrame({'released_state':observed_states,'alpha_s':alpha,'mu':coef[0]}),patient_effects=pd.DataFrame({'patient':patients,'gamma_i':gamma}),residuals=residuals,design=mat)
