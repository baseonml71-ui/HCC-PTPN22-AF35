from common import *
from itertools import combinations
from math import comb
feas=read(OUT/(P+'_FEASIBILITY.tsv'));data=read(OUT/'inputs/PATIENT_COMPONENTS_AFTER_RESPONSE_JOIN.tsv');grid=read(OUT/'AF35_PATIENT_STATE_PROPORTIONS.tsv');state=read(OUT/'AF35_STATE_EFFECTS_RESPONSE_BLIND.tsv');cap=read(OUT/'AF35_STATE_OCCUPANCY_CELLCOUNT_SENSITIVITY_DESIGN.tsv')
keys=['M','O','W','X'];cols=[k+'_i' for k in keys]
def g_effect(a,b):
 n=a.shape[-2]+b.shape[-2];v=((a.shape[-2]-1)*a.var(axis=-2,ddof=1)+(b.shape[-2]-1)*b.var(axis=-2,ddof=1))/(n-2)
 valid=(v>0)&~((np.ptp(a,axis=-2)==0)&(np.ptp(b,axis=-2)==0))
 return np.divide((1-3/(4*(n-2)-1))*(a.mean(axis=-2)-b.mean(axis=-2)),np.sqrt(v),out=np.full_like(v,np.nan),where=valid)
effects=[];boots=[];summ=[];hedges=[];perm=[];perm_detail=[];contributions=[];decisions=[];sensitivity=[];identity=[]
for compartment in feas[feas.status=='PASS'].compartment:
 d=data[(data.compartment==compartment)&data.response.notna()].sort_values('patient');patients=d.patient.tolist();x=d[cols].to_numpy();r=d.response.to_numpy()==1;n=len(x);nr=int(r.sum());nn=n-nr
 delta=x[r].mean(axis=0)-x[~r].mean(axis=0);diff=delta[1]-delta[2];g=g_effect(x[r],x[~r])
 rng=np.random.default_rng(2026090810);ri=np.flatnonzero(r);ni=np.flatnonzero(~r)
 bix=np.column_stack([rng.choice(ri,size=(10000,nr),replace=True),rng.choice(ni,size=(10000,nn),replace=True)])
 save(pd.DataFrame(bix,columns=[f'slot_{j+1}' for j in range(n)]),f'inputs/{compartment}_BOOTSTRAP_INDICES_ZERO_BASED.tsv');save({'position_zero_based':range(n),'patient':patients,'response':r.astype(int)},f'inputs/{compartment}_BOOTSTRAP_PATIENT_ORDER.tsv')
 a=x[bix[:,:nr]];b=x[bix[:,nr:]];bs=a.mean(axis=1)-b.mean(axis=1);bd=bs[:,1]-bs[:,2];bg=g_effect(a,b)
 boot=pd.DataFrame(bs,columns=['Delta_'+k for k in keys]);boot['Delta_O_minus_W']=bd;boot.insert(0,'draw',np.arange(1,10001));boot.insert(0,'compartment',compartment);boots.append(boot)
 raw_error=float(abs(delta[0]-delta[1:].sum()));boot_error=float(np.max(abs(bs[:,0]-bs[:,1:].sum(axis=1))));assert max(raw_error,boot_error)<1e-12
 for j,key in enumerate(keys+['O_minus_W']):
  vv=bs[:,j] if j<4 else bd;value=delta[j] if j<4 else diff;lo,hi=np.quantile(vv,[.025,.975])
  effects.append(dict(compartment=compartment,component=key,R_n=nr,NR_n=nn,mean_R=float(x[r,j].mean()) if j<4 else float((x[r,1]-x[r,2]).mean()),mean_NR=float(x[~r,j].mean()) if j<4 else float((x[~r,1]-x[~r,2]).mean()),raw_R_minus_NR=value,bootstrap_ci_low=lo,bootstrap_ci_high=hi,units='FROZEN_AF35_SCORE'))
  summ.append(dict(compartment=compartment,component=key,draws=10000,observed=value,bootstrap_median=float(np.median(vv)),ci_low=lo,ci_high=hi,interval='PATIENT_GROUP_STRATIFIED_PERCENTILE_95',frozen_component_model_not_refit=True))
  if j<4:
   good=bg[:,j][np.isfinite(bg[:,j])];glo,ghi=np.quantile(good,[.025,.975]) if len(good) else [np.nan,np.nan]
   hedges.append(dict(compartment=compartment,component=key,R_n=nr,NR_n=nn,hedges_g=g[j],bootstrap_ci_low=glo,bootstrap_ci_high=ghi,finite_bootstrap_draws=len(good),nonfinite_draws=10000-len(good),interpretation='SECONDARY_NON_ADDITIVE_STANDARDIZED_EFFECT'))
 lo,hi=np.quantile(bd,[.025,.975]);pattern='OCCUPANCY_LEADING' if diff>0 and lo>0 else ('WITHIN_STATE_LEADING' if diff<0 and hi<0 else 'MIXED_OR_IMPRECISE')
 interaction=bool(abs(delta[3])>max(abs(delta[1]),abs(delta[2])));bridge='POSITIVE_DIRECTION' if delta[0]>0 else 'DIRECTION_NOT_REPLICATED';final=bridge if delta[0]<=0 else ('INTERACTION_DOMINANT' if interaction else pattern)
 decisions.append(dict(compartment=compartment,occupancy_vs_within_state=pattern,interaction_dominant=interaction,response_bridge=bridge,final_scientific_pattern=final,Delta_O_minus_W=diff,comparison_ci_low=lo,comparison_ci_high=hi))
 pi=grid[grid.compartment==compartment].pivot(index='patient',columns='released_state',values='p_is').loc[patients];st=state[state.compartment==compartment].sort_values('state_order')
 for row in st.itertuples():
  v=pi[row.released_state].to_numpy();dp=float(v[r].mean()-v[~r].mean());cs=dp*row.alpha_s if row.estimated else 0.
  contributions.append(dict(compartment=compartment,state_order=row.state_order,released_state=row.released_state,alpha_s=row.alpha_s,estimated=row.estimated,mean_p_R=float(v[r].mean()),mean_p_NR=float(v[~r].mean()),Delta_p=dp,C_s=cs))
 cerror=abs(sum(v['C_s'] for v in contributions if v['compartment']==compartment)-delta[1]);assert cerror<1e-12
 identity.append(dict(compartment=compartment,patient_reconstruction_max_error=float(d.reconstruction_error.abs().max()),group_reconstruction_error=raw_error,bootstrap_reconstruction_max_error=boot_error,state_contribution_sum_error=cerror,status='PASS'))
 assignments=list(combinations(range(n),nr));pdeltas=[]
 for number,assignment in enumerate(assignments):
  rr=np.zeros(n,dtype=bool);rr[list(assignment)]=True;dv=x[rr].mean(axis=0)-x[~rr].mean(axis=0);pdeltas.append(dv)
  perm_detail.append(dict(compartment=compartment,assignment_id=number+1,responder_patients=';'.join(patients[i] for i in assignment),Delta_M=dv[0],Delta_O=dv[1],Delta_W=dv[2]))
 pdeltas=np.array(pdeltas);assert len(assignments)==comb(n,nr)
 for j,key in enumerate(keys[:3]):
  extreme=int((np.abs(pdeltas[:,j])>=abs(delta[j])-1e-12).sum());perm.append(dict(compartment=compartment,component=key,observed=delta[j],R_n=nr,NR_n=nn,unique_assignments=len(assignments),extreme_assignments=extreme,exact_two_sided_P=extreme/len(assignments),tail_rule='ABS_PERM_GE_ABS_OBS_TOL_1e-12',supportive_only=True))
 cd=cap[cap.compartment==compartment].iloc[0];sensitivity.append(dict(compartment=compartment,status=cd.status,cap=int(cd.cap),retained_cells=int(cd.retained_cells),original_cells=int(cd.original_cells),retained_fraction=cd.retained_fraction,reason='超过结局前预定90%细胞丢失界限；不改cap、不重抽救援'))
 if cd.status=='RUN_FROZEN_BEFORE_RESPONSE':
  cx=read(OUT/'inputs/CAPPED_PATIENT_COMPONENTS_RESPONSE_BLIND.tsv').query('compartment == @compartment').set_index('patient').loc[patients,cols].to_numpy();de=cx[r].mean(axis=0)-cx[~r].mean(axis=0);bv=cx[bix[:,:nr]].mean(axis=1)-cx[bix[:,nr:]].mean(axis=1);ll,hh=np.quantile(bv[:,1]-bv[:,2],[.025,.975]);sensitivity[-1].update({f'Delta_{k}':de[j] for j,k in enumerate(keys)});sensitivity[-1].update(O_minus_W_ci_low=ll,O_minus_W_ci_high=hh)
for n,rows in [(P+'_CLINICAL_DECOMPOSITION.tsv',effects),(P+'_BOOTSTRAP_SUMMARY.tsv',summ),(P+'_HEDGES_G.tsv',hedges),(P+'_EXACT_PERMUTATION.tsv',perm),(P+'_STATE_CONTRIBUTIONS.tsv',contributions),(P+'_SCIENTIFIC_INTERPRETATION.tsv',decisions),(P+'_CELLCOUNT_SENSITIVITY.tsv',sensitivity),('QA/EXACT_IDENTITY_VERIFICATION.tsv',identity),('QA/EXACT_PERMUTATION_ASSIGNMENTS.tsv',perm_detail)]:save(rows,n)
save(pd.concat(boots) if boots else pd.DataFrame(columns=['compartment','draw','Delta_M','Delta_O','Delta_W','Delta_X','Delta_O_minus_W']),P+'_BOOTSTRAP.tsv')
js('QA/CLINICAL_COMPLETION_LOCK.json',dict(utc=now(),passing_compartments=feas[feas.status=='PASS'].compartment.tolist(),scientific_patterns=decisions,hashes={p.name:sha(p) for p in OUT.glob(P+'_*.tsv')}))
print(pd.DataFrame(effects).to_string(index=False));print(pd.DataFrame(decisions).to_string(index=False));print(pd.DataFrame(hedges).to_string(index=False));print(pd.DataFrame(perm).to_string(index=False))
