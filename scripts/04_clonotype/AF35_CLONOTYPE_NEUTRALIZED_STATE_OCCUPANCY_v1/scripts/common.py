from pathlib import Path
import datetime,hashlib,json,os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
ROOT=Path('.');OUT=Path(__file__).resolve().parents[1]
K=ROOT/'AF35_CLINICAL_STATE_OCCUPANCY_DECOMPOSITION_v1';D=ROOT/'GSE235863_CLONE_STATE_DECOMPOSITION_v1';P='AF35_CLONE_NEUTRALIZATION'
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
def neutralize(cells,identity,alpha,states):
 d=cells.rename(columns={identity:'clone_identity'}).copy()
 counts=d.groupby(['patient','clone_identity','released_state'],sort=True).size().rename('n_ics').reset_index()
 counts['n_ic']=counts.groupby(['patient','clone_identity']).n_ics.transform('sum');counts['f_ics']=counts.n_ics/counts.n_ic;counts['alpha_s']=counts.released_state.map(alpha);counts['f_alpha']=counts.f_ics*counts.alpha_s
 potentials=counts.groupby(['patient','clone_identity'],sort=True).agg(n_ic=('n_ic','first'),A_ic=('f_alpha','sum'),occupied_states=('released_state','nunique'),sum_f=('f_ics','sum')).reset_index()
 assert np.max(abs(potentials.sum_f-1))<1e-14
 proportions=[];components=[];coupling=[];concentration=[]
 for pat,g in potentials.groupby('patient',sort=True):
  cc=len(g);n=int(g.n_ic.sum());r=counts[counts.patient==pat];pr=r.groupby('released_state').agg(n_cells=('n_ics','sum'),sum_clone_fraction=('f_ics','sum')).reindex(states).fillna(0)
  pr['p_cell']=pr.n_cells/n;pr['p_clone']=pr.sum_clone_fraction/cc;pr['alpha_s']=[alpha[s] for s in states];pr['patient']=pat;pr['state_order']=range(1,21);proportions.append(pr.rename_axis('released_state').reset_index())
  oc=float(np.dot(pr.p_cell,pr.alpha_s));oe=float(np.dot(pr.p_clone,pr.alpha_s));ds=oc-oe
  oc2=float(np.dot(g.n_ic,g.A_ic)/n);oe2=float(g.A_ic.mean());cov=float(np.mean((g.n_ic-g.n_ic.mean())*(g.A_ic-g.A_ic.mean())));ratio=cov/g.n_ic.mean()
  di=(g.n_ic/n-1/cc)*g.A_ic;potentials.loc[g.index,'d_ic']=di
  components.append(dict(patient=pat,resolved_cells=n,clone_count=cc,O_CELL=oc,O_CLONE=oe,D_SIZE=ds,O_CELL_clone_formula=oc2,O_CLONE_clone_formula=oe2,population_cov_n_A=cov,mean_local_clone_size=g.n_ic.mean(),covariance_over_mean=ratio,cell_formula_error=oc-oc2,clone_formula_error=oe-oe2,patient_identity_error=oc-(oe+ds),covariance_identity_error=ds-ratio,clone_contribution_sum_error=ds-di.sum()))
  finite=g.n_ic.nunique()>1 and g.A_ic.nunique()>1
  coupling.append(dict(patient=pat,clone_count=cc,median_local_clone_size=float(g.n_ic.median()),maximum_local_clone_size=int(g.n_ic.max()),spearman_log1p_size_potential=float(spearmanr(np.log1p(g.n_ic),g.A_ic).statistic) if finite else np.nan,status='DEFINED' if finite else 'UNDEFINED_CONSTANT_VECTOR'))
  absolute=np.abs(di);mass=float(absolute.sum());pos=di[di>0];neg=di[di<0]
  concentration.append(dict(patient=pat,clone_count=cc,largest_positive_d=float(pos.max()) if len(pos) else 0.,largest_negative_d=float(neg.min()) if len(neg) else 0.,total_absolute_d_mass=mass,top3_absolute_mass_fraction=float(np.sort(absolute)[-3:].sum()/mass) if mass else np.nan,sum_d=float(di.sum()),D_SIZE=ds))
 comp=pd.DataFrame(components)
 assert comp[['cell_formula_error','clone_formula_error','patient_identity_error','covariance_identity_error','clone_contribution_sum_error']].abs().to_numpy().max()<1e-12
 return counts,potentials,pd.concat(proportions,ignore_index=True),comp,pd.DataFrame(coupling),pd.DataFrame(concentration)
