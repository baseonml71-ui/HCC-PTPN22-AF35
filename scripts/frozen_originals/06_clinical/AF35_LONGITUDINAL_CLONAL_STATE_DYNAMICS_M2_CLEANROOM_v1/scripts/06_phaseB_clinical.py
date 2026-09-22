from common import *
import itertools, shutil
assert not any(B.iterdir()) and not any(C.iterdir())
lock=verify_phasea();assert sha(OUT/(P+'_PHASEA_COMPLETE_SHA256_MANIFEST.tsv'))==lock['manifest_sha256']
lockhash=sha(OUT/(P+'_PHASEA_FREEZE_LOCK.json'))
js(Q/'PHASEB_START_EVENT.json',dict(utc=now(),phaseA_lock_sha256=lockhash,phaseA_freeze_utc=lock['utc'],phaseA_hash_verification='PASS',first_clinical_read_follows_this_event=True))
source=K/'inputs/RESPONSE_LABEL_MAP_AFTER_FREEZE.tsv'
manifest=read(K/'COMPLETE_SHA256_MANIFEST.tsv');e=manifest[manifest.path.str.replace('\\','/',regex=False).eq('inputs/RESPONSE_LABEL_MAP_AFTER_FREEZE.tsv')].sha256.iloc[0]
assert sha(source)==e
lab=read(source,usecols=['patient','response_group']).rename(columns={'response_group':'response'})
assert lab.patient.is_unique and set(lab.response)=={'R','NR'}
save(lab,B/'RESPONSE.tsv')
delta=read(F/(P+'_PATIENT_LONGITUDINAL_CHANGES_RESPONSE_BLIND.tsv'));assert delta.patient.is_unique
joined=delta.merge(lab,on='patient',how='left',validate='one_to_one');assert joined.response.notna().all() and set(joined.patient)==set(delta.patient)
for c in delta.columns:assert np.array_equal(joined[c],delta[c]),'Frozen molecular metric changed'
save(joined,C/'PATIENT_CHANGES_AFTER_JOIN.tsv')
save(lab.assign(included=lab.patient.isin(delta.patient)),C/'RESPONSE_JOIN_COVERAGE.tsv')
counts=joined.response.value_counts();assert min(counts.get('R',0),counts.get('NR',0))>=2,'NOT_EVALUABLE_LOW_PAIRED_COVERAGE'
write(OUT/(P+'_PHASEB_RESPONSE_SOURCE_AUDIT.md'),f'''# 阶段 B 响应来源

唯一患者响应来源：`{source}`。SHA-256：`{e}`，与 K 完整清单一致。仅提取 patient、response_group，并将字段名改为 response，保留既有 R/NR 值，未载入其他临床协变量。

阶段 A 完成于 {lock['utc']}，接入前全部冻结哈希通过；首次临床读取之前的阶段 B 开始事件见 QA/PHASEB_START_EVENT.json。分子患者 {len(joined)} 人，R={counts['R']}，NR={counts['NR']}，无重复、无未匹配、无额外临床患者加入。源中仅 post 的额外患者未进入分子配对。

L 的群体基线数值作为另一个已冻结的组水平参照，仅在本阶段读取，不作为患者响应映射，不重算。
''')
columns=['CHANGE_O_CELL','CHANGE_O_CLONE','CHANGE_D_SIZE'];names=['THETA_CELL','THETA_CLONE','THETA_SIZE']
def infer(table,cols,names0,seed,prefix):
 v=table[cols].to_numpy();ir=np.flatnonzero(table.response.eq('R'));inn=np.flatnonzero(table.response.eq('NR'));rng=np.random.default_rng(seed)
 ri=rng.choice(ir,size=(10000,len(ir)),replace=True);ni=rng.choice(inn,size=(10000,len(inn)),replace=True)
 draws=v[ri].mean(axis=1)-v[ni].mean(axis=1);est=v[ir].mean(axis=0)-v[inn].mean(axis=0)
 use=list(names0)
 if len(cols)==3:draws=np.column_stack([draws,draws[:,1]-draws[:,2]]);est=np.append(est,est[1]-est[2]);use+=['THETA_BALANCE']
 save(dict(patient_index=np.arange(len(table)),patient=table.patient,response=table.response),C/(prefix+'_PATIENT_ORDER.tsv'))
 save(pd.DataFrame(np.column_stack([ri,ni]),columns=[f'R_{i}' for i in range(len(ir))]+[f'NR_{i}' for i in range(len(inn))]),C/(prefix+'_BOOTSTRAP_INDICES.tsv'))
 boot=pd.DataFrame(draws,columns=use);boot.insert(0,'draw',np.arange(10000));save(boot,C/(prefix+'_BOOTSTRAP.tsv'))
 distributions=[]
 for ri0 in itertools.combinations(range(len(table)),len(ir)):
  mask=np.zeros(len(table),dtype=bool);mask[list(ri0)]=True;ev=v[mask].mean(axis=0)-v[~mask].mean(axis=0)
  distributions.append(dict(R_indices=','.join(map(str,ri0)),**dict(zip(names0,ev))))
 perms=pd.DataFrame(distributions);save(perms,C/(prefix+'_EXACT_PERMUTATION_DISTRIBUTION.tsv'))
 results=[]
 for j,n in enumerate(use):
  lo,hi=np.quantile(draws[:,j],[.025,.975]);row=dict(component=n,estimate=est[j],bootstrap_ci_low=lo,bootstrap_ci_high=hi,bootstrap_draws=10000,seed=seed,R_n=len(ir),NR_n=len(inn))
  row['exact_two_sided_p']=float(np.mean(abs(perms[n].to_numpy())>=abs(est[j])-1e-12)) if n in perms else np.nan
  row['exact_assignments']=len(perms) if n in perms else np.nan;results.append(row)
 result=pd.DataFrame(results);save(result,C/(prefix+'_BOOTSTRAP_SUMMARY.tsv'))
 save(result[result.component.isin(names0)][['component','estimate','exact_two_sided_p','exact_assignments']],C/(prefix+'_EXACT_PERMUTATION.tsv'))
 return result,boot
res,boot=infer(joined,columns,names,2026090821,P)
assert abs(res.estimate.iloc[0]-res.estimate.iloc[1]-res.estimate.iloc[2])<1e-12
assert abs(boot.THETA_CELL-boot.THETA_CLONE-boot.THETA_SIZE).max()<1e-12
contrast=res.copy()
for j,n in enumerate(names):
 contrast.loc[j,'mean_R_change']=joined.loc[joined.response.eq('R'),columns[j]].mean();contrast.loc[j,'mean_NR_change']=joined.loc[joined.response.eq('NR'),columns[j]].mean()
save(contrast,C/(P+'_CLINICAL_LONGITUDINAL_CONTRASTS.tsv'))
baseline_source=L/'AF35_CLONE_NEUTRALIZATION_CLINICAL_DECOMPOSITION.tsv'
lm=read(L/'COMPLETE_SHA256_MANIFEST.tsv');assert sha(baseline_source)==lm[lm.path.eq(baseline_source.name)].sha256.iloc[0]
base=read(baseline_source);base=base[base.component.isin(['O_CELL','O_CLONE','D_SIZE'])].copy();save(base,B/'L_FROZEN_GROUP_BASELINE.tsv')
baseline=float(base[base.component.eq('O_CELL')].raw_R_minus_NR.iloc[0]);theta=float(res.estimate.iloc[0])
pattern='PRE_EXISTING_POSITIVE_WITH_REINFORCEMENT' if baseline>0 and theta>0 else 'PRE_EXISTING_POSITIVE_WITHOUT_REINFORCEMENT' if baseline>0 else 'LONGITUDINAL_POSITIVE_WITH_WEAK_BASELINE' if theta>0 else 'NO_CONSISTENT_POSITIVE_LONGITUDINAL_PATTERN'
bal=res[res.component.eq('THETA_BALANCE')].iloc[0];dominance='CLONE_EQUAL_CHANGE_LEADING' if bal.bootstrap_ci_low>0 else 'CLONE_SIZE_CHANGE_LEADING' if bal.bootstrap_ci_high<0 else 'MIXED_OR_IMPRECISE'
time=read(F/(P+'_PATIENT_TIMEPOINT_COMPONENTS_RESPONSE_BLIND.tsv')).merge(lab,on='patient',validate='many_to_one')
des=[]
for group in ['R','NR']:
 for col,chg in zip(['O_CELL','O_CLONE','D_SIZE'],columns):
  a=time[time.response.eq(group)&time.timepoint.eq('pre')][col];b=time[time.response.eq(group)&time.timepoint.eq('post')][col];v=joined[joined.response.eq(group)][chg]
  des.append(dict(response=group,component=col,n=len(v),pre_mean=a.mean(),later_mean=b.mean(),change_mean=v.mean(),change_median=v.median(),positive_change_patients=int((v>0).sum())))
save(des,C/'WITHIN_GROUP_LONGITUDINAL_DESCRIPTION.tsv')
sc=read(F/(P+'_STATE_CONTRIBUTIONS_RESPONSE_BLIND.tsv')).merge(lab,on='patient',validate='many_to_one');sr=[]
for group in ['ALL','R','NR','R_MINUS_NR']:
 for (order,state),g in sc.groupby(['state_order','released_state'],sort=True):
  row=dict(group=group,state_order=order,released_state=state)
  for col in ['CELL_CONTRIBUTION','CLONE_CONTRIBUTION']:
   row[col]=g[col].mean() if group=='ALL' else g[g.response.eq('R')][col].mean()-g[g.response.eq('NR')][col].mean() if group=='R_MINUS_NR' else g[g.response.eq(group)][col].mean()
  sr.append(row)
state=pd.DataFrame(sr);save(state,C/'STATE_CONTRIBUTIONS.tsv')
assert abs(state[state.group.eq('R_MINUS_NR')].CELL_CONTRIBUTION.sum()-theta)<1e-12
assert abs(state[state.group.eq('R_MINUS_NR')].CLONE_CONTRIBUTION.sum()-res.estimate.iloc[1])<1e-12
shared=read(F/(P+'_SHARED_CLONE_METRICS_RESPONSE_BLIND.tsv'));shared=shared[shared.shift_eligible].merge(lab,on='patient',validate='one_to_one');sg=shared.response.value_counts()
sharedgate='NOT_EVALUABLE_LOW_SHARED_CLONE_COVERAGE'
if len(shared)>=5 and min(sg.get('R',0),sg.get('NR',0))>=2:
 sharedres,_=infer(shared,['SHARED_CLONE_MEAN_SHIFT'],['THETA_SHARED'],2026090823,P+'_SHARED')
 save(sharedres,C/(P+'_SHARED_CLONE_CLINICAL.tsv'));sharedgate='PASS_SECONDARY_ONLY'
elif len(shared)>=5:sharedgate='NOT_EVALUABLE_LOW_SHARED_RESPONSE_GROUP_COVERAGE'
verify_phasea();assert sha(OUT/(P+'_PHASEA_FREEZE_LOCK.json'))==lockhash
js(C/'CLINICAL_COMPLETION.json',dict(utc=now(),patients=len(joined),R=int(counts['R']),NR=int(counts['NR']),pattern=pattern,balance=dominance,shared_gate=sharedgate,phaseA_post_join_hashes='PASS',baseline_source=str(baseline_source),baseline_source_sha256=sha(baseline_source),baseline_recomputed=False,baseline_same_molecular_patient_universe=set(read(A/'L_IDENTITY_REFERENCE.tsv').patient)==set(joined.patient),maximum_theta_identity_error=float(abs(res.estimate.iloc[0]-res.estimate.iloc[1]-res.estimate.iloc[2])),maximum_bootstrap_identity_error=float(abs(boot.THETA_CELL-boot.THETA_CLONE-boot.THETA_SIZE).max())))
print(contrast.to_string(index=False),flush=True);print(pattern,dominance,sharedgate,flush=True)
