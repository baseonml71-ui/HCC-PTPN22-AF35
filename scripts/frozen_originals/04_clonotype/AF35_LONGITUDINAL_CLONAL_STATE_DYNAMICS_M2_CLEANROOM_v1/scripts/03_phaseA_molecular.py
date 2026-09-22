from common import *
assert not (F/(P+'_PATIENT_TIMEPOINT_COMPONENTS_RESPONSE_BLIND.tsv')).exists()
events=install_guard()
# Attempted opens are rejected before the operating system opens the files.
for forbidden in [ROOT/'AF35_LONGITUDINAL_CLONAL_STATE_DYNAMICS_v1/README.md',L/'AF35_CLONE_NEUTRALIZATION_CLINICAL_DECOMPOSITION.tsv',B/'RESPONSE.tsv']:
 try:
  open(forbidden,'rb');raise AssertionError('Clean room access control failed')
 except PermissionError:pass
x=read(A/'MOLECULAR_CELLS.tsv'); alpha=read(A/'K_FROZEN_ALPHA.tsv'); states=read(A/'RELEASED_STATE_ORDER.tsv').released_state.tolist()
assert alpha.released_state.tolist()==states and len(states)==20 and set(alpha.compartment)=={'P'}
av=alpha.alpha_s.to_numpy(); amap=dict(zip(states,av));assert set(x.released_state).issubset(states)
x=x[x.compartment.eq('P')].copy()
sample=x.groupby(['patient','sample','timepoint'],sort=True).agg(cell_n=('barcode','size'),clone_n=('clonotype_key','nunique'),state_n=('released_state','nunique'),h5_sample=('h5_sample','first'),h5_sample_n=('h5_sample','nunique')).reset_index()
patients=[]; sample['pairing_status']=''
for pat,g in sample.groupby('patient',sort=True):
 counts=g.groupby('timepoint')['sample'].nunique()
 unique=counts.get('pre',0)==1 and counts.get('post',0)==1
 native=(g['sample']==g.h5_sample).all() and g.h5_sample_n.eq(1).all()
 status='UNIQUELY_RESOLVED_OUTCOME_FREE' if unique and native else 'PAIRING_AMBIGUOUS_OUTCOME_FREE' if counts.get('pre',0)>1 or counts.get('post',0)>1 or not native else 'EXCLUDED_MISSING_PRE_OR_POST'
 sample.loc[g.index,'pairing_status']=status
 if status=='UNIQUELY_RESOLVED_OUTCOME_FREE':patients.append(pat)
save(sample,F/(P+'_OUTCOME_FREE_SAMPLE_MANIFEST.tsv'))
assert len(patients)>=6, 'AF35_M2_MOLECULAR_COVERAGE_GATE=NOT_EVALUABLE_LOW_PAIRED_COVERAGE'
x=x[x.patient.isin(patients)]
assert (x.clonotype_key==x.patient+'::'+x['clone.id'].astype(str)).all()
assert x.barcode.is_unique and x.AF35.notna().all()
save(x,F/'PAIRED_MOLECULAR_CELLS.tsv')
ct=x.groupby(['patient','timepoint','clonotype_key','released_state'],sort=True).size().rename('n_icst').reset_index()
ct['n_ict']=ct.groupby(['patient','timepoint','clonotype_key']).n_icst.transform('sum')
ct['f_icst']=ct.n_icst/ct.n_ict;ct['alpha_s']=ct.released_state.map(amap);ct['f_alpha']=ct.f_icst*ct.alpha_s
pot=ct.groupby(['patient','timepoint','clonotype_key'],sort=True).agg(n_ict=('n_ict','first'),A_ict=('f_alpha','sum'),sum_f=('f_icst','sum')).reset_index()
save(ct,F/'CLONE_STATE_DISTRIBUTIONS.tsv');save(pot,F/'CLONE_STATE_POTENTIAL.tsv')
rows=[];props=[]
for (pat,t),g in pot.groupby(['patient','timepoint'],sort=True):
 r=ct[ct.patient.eq(pat)&ct.timepoint.eq(t)]
 pr=r.groupby('released_state').agg(n_cells=('n_icst','sum'),sum_f=('f_icst','sum')).reindex(states).fillna(0)
 n=int(g.n_ict.sum());cc=len(g);pr['p_cell']=pr.n_cells/n;pr['p_clone']=pr.sum_f/cc
 oc=float(pr.p_cell.to_numpy()@av);oe=float(pr.p_clone.to_numpy()@av);ds=oc-oe
 oc2=float(g.n_ict.to_numpy()@g.A_ict.to_numpy()/n);oe2=float(g.A_ict.mean());cov=float(np.mean((g.n_ict-g.n_ict.mean())*(g.A_ict-g.A_ict.mean()))/g.n_ict.mean())
 rows.append(dict(patient=pat,timepoint=t,compartment='P',cell_n=n,clone_n=cc,state_n=int((pr.n_cells>0).sum()),O_CELL=oc,O_CLONE=oe,D_SIZE=ds,identity_error=oc-oe-ds,cell_clone_formula_error=oc-oc2,clone_formula_error=oe-oe2,covariance_identity_error=ds-cov))
 pr['patient']=pat;pr['timepoint']=t;pr['state_order']=range(1,21);pr['alpha_s']=av;props.append(pr.rename_axis('released_state').reset_index())
comp=pd.DataFrame(rows);pr=pd.concat(props,ignore_index=True)
assert comp.filter(like='error').abs().to_numpy().max()<1e-12
save(comp,F/(P+'_PATIENT_TIMEPOINT_COMPONENTS_RESPONSE_BLIND.tsv'));save(pr,F/'PATIENT_STATE_PROPORTIONS.tsv')
deltas=[];statechanges=[]
for pat in patients:
 c=comp[comp.patient.eq(pat)].set_index('timepoint');row=dict(patient=pat)
 for n in ['O_CELL','O_CLONE','D_SIZE']:row['CHANGE_'+n]=c.loc['post',n]-c.loc['pre',n]
 row['identity_error']=row['CHANGE_O_CELL']-row['CHANGE_O_CLONE']-row['CHANGE_D_SIZE'];deltas.append(row)
 q=pr[pr.patient.eq(pat)];pre=q[q.timepoint.eq('pre')].set_index('released_state').loc[states];post=q[q.timepoint.eq('post')].set_index('released_state').loc[states]
 for j,s in enumerate(states):statechanges.append(dict(patient=pat,state_order=j+1,released_state=s,CELL_CONTRIBUTION=(post.loc[s,'p_cell']-pre.loc[s,'p_cell'])*av[j],CLONE_CONTRIBUTION=(post.loc[s,'p_clone']-pre.loc[s,'p_clone'])*av[j]))
delta=pd.DataFrame(deltas); sc=pd.DataFrame(statechanges)
assert delta.identity_error.abs().max()<1e-12
save(delta,F/(P+'_PATIENT_LONGITUDINAL_CHANGES_RESPONSE_BLIND.tsv'));save(sc,F/(P+'_STATE_CONTRIBUTIONS_RESPONSE_BLIND.tsv'))
shared=[];coverage=[];metrics=[];memory=[]
for pi,pat in enumerate(patients):
 g=pot[pot.patient.eq(pat)]
 pre=g[g.timepoint.eq('pre')].set_index('clonotype_key');post=g[g.timepoint.eq('post')].set_index('clonotype_key')
 ids=sorted(set(pre.index)&set(post.index));n=len(ids);pa=pre.loc[ids];pb=post.loc[ids]
 for key in ids:shared.append(dict(patient=pat,clonotype_key=key,A_pre=pa.loc[key,'A_ict'],A_later=pb.loc[key,'A_ict'],CHANGE_A=pb.loc[key,'A_ict']-pa.loc[key,'A_ict'],pre_n=int(pa.loc[key,'n_ict']),later_n=int(pb.loc[key,'n_ict'])))
 coverage.append(dict(patient=pat,shared_clone_n=n,pre_clone_fraction=n/len(pre),later_clone_fraction=n/len(post),pre_shared_cells=int(pa.n_ict.sum()),later_shared_cells=int(pb.n_ict.sum()),pre_cell_representation=pa.n_ict.sum()/pre.n_ict.sum(),later_cell_representation=pb.n_ict.sum()/post.n_ict.sum(),shift_eligible=n>=3,memory_eligible=n>=5))
 shift=float(np.mean(pb.A_ict.to_numpy()-pa.A_ict.to_numpy())) if n>=3 else np.nan
 metrics.append(dict(patient=pat,shared_clone_n=n,SHARED_CLONE_MEAN_SHIFT=shift,shift_eligible=n>=3))
 if n>=5:
  aa=rankdata(np.round(pa.A_ict.to_numpy(),12));bb=rankdata(np.round(pb.A_ict.to_numpy(),12));aa-=aa.mean();bb-=bb.mean();den=np.linalg.norm(aa)*np.linalg.norm(bb)
  seed=2026090822+pi;rng=np.random.default_rng(seed);indices=np.stack([rng.permutation(n) for _ in range(1000)])
  save(pd.DataFrame(indices),F/f'MEMORY_INDICES_{pat}.tsv')
  if den>0:
   obs=float(aa@bb/den);dist=(bb[indices]@aa)/den
   mid=float((np.sum(dist<obs-1e-12)+.5*np.sum(abs(dist-obs)<=1e-12))/1000)
   memory.append(dict(patient=pat,shared_clone_n=n,SAME_CLONE_STATE_POTENTIAL_CORRELATION=obs,reference_q025=np.quantile(dist,.025),reference_median=np.median(dist),reference_q975=np.quantile(dist,.975),reference_fraction_below=float(np.mean(dist<obs-1e-12)),reference_fraction_equal=float(np.mean(abs(dist-obs)<=1e-12)),reference_midrank_position=mid,permutations=1000,seed=seed,status='DEFINED'))
   save(dict(draw=np.arange(1000),correlation=dist),F/f'MEMORY_REFERENCE_{pat}.tsv')
  else:memory.append(dict(patient=pat,shared_clone_n=n,SAME_CLONE_STATE_POTENTIAL_CORRELATION=np.nan,permutations=1000,seed=seed,status='UNDEFINED_CONSTANT_VECTOR'))
save(shared,F/'SHARED_CLONE_POTENTIAL_VECTORS.tsv');save(coverage,F/'SHARED_CLONE_COVERAGE.tsv');save(metrics,F/(P+'_SHARED_CLONE_METRICS_RESPONSE_BLIND.tsv'));save(memory,F/'SAME_CLONE_MEMORY_RESPONSE_BLIND.tsv')
js(F/'MOLECULAR_COMPLETION.json',dict(utc=now(),paired_patients=patients,paired_n=len(patients),P27_status=sample[sample.patient.eq('P27')].pairing_status.unique().tolist(),shared_eligible_patients=int(pd.DataFrame(metrics).shift_eligible.sum()),source_native_P_only=True,global_clone_size_read=False,clinical_data_accessed=False,alpha_refit=False,maximum_identity_error=float(max(comp.filter(like='error').abs().to_numpy().max(),delta.identity_error.abs().max()))))
save(events,F/'PHASEA_FILE_ACCESS_LOG.tsv')
print('Phase A molecular complete:',len(patients),'pairs; no outcome data accessed',flush=True)
