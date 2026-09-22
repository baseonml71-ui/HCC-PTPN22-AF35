import importlib,json
import numpy as np
import pandas as pd
f=importlib.import_module('00_initialize');OUT=f.OUT;P=f.PREFIX
lock=json.loads((OUT/(P+'_FREEZE_LOCK.json')).read_text(encoding='utf-8'))
for n,h in lock['hashes'].items():assert f.sha(OUT/n)==h
for n,h in lock['source_hashes'].items():assert f.sha(n)==h
assert lock['outcome_labels_loaded'] is False and lock['all_LOO_residuals_frozen'] is True
assert not (OUT/(P+'_OUTCOME_JOIN_TIMESTAMP.txt')).exists()
f.write(P+'_OUTCOME_JOIN_TIMESTAMP.txt',f.now()+'\n')
def hg(x,y):
 r=x[y==1];nr=x[y==0];n=len(x)
 sd=np.sqrt(((len(r)-1)*r.var(ddof=1)+(len(nr)-1)*nr.var(ddof=1))/(n-2))
 return (1-3/(4*(n-2)-1))*(r.mean()-nr.mean())/sd
def bhg(x,idx,nr):
 v=x[idx];r=v[:,:nr];n=v[:,nr:];N=v.shape[1]
 sd=np.sqrt(((r.shape[1]-1)*r.var(axis=1,ddof=1)+(n.shape[1]-1)*n.var(axis=1,ddof=1))/(N-2))
 return (1-3/(4*(N-2)-1))*(r.mean(axis=1)-n.mean(axis=1))/sd
values=f.read(OUT/(P+'_RESPONSE_BLIND_PATIENT_VALUES.tsv'))
folds=f.read(OUT/(P+'_LOO_RESPONSE_BLIND_PATIENT_VALUES.tsv'))
effects=[];boots=[];summaries=[];loos=[];looss=[];bindings=[];joined=[]
for k,c in enumerate(f.COHORTS):
 path=f.A/f'{c}_ANALYSIS_PATIENT_VALUES.tsv';target=f.A/'AF35_FULL_SCORE_EFFECT_REPRODUCTION.tsv'
 bindings.extend([dict(path=str(p),sha256=f.sha(p),first_content_read_after_freeze=True) for p in [path,target]])
 d=values[values.cohort==c].copy();old=f.read(path,usecols=['patient_id','response','AF35'])
 assert set(d.patient_id)==set(old.patient_id) and old.patient_id.is_unique
 old=old.set_index('patient_id').loc[d.patient_id];assert np.max(abs(old.AF35.to_numpy()-d.AF35.to_numpy()))<1e-12
 y=old.response.to_numpy(dtype=int);assert set(y)=={0,1};d['response']=y;joined.append(d)
 x=d.AF35.to_numpy();res=d.AF35_residual.to_numpy();a=hg(x,y);b=hg(res,y)
 expected=f.read(target).query('cohort==@c').frozen_hedges_g.item();assert abs(a-expected)<1e-10
 effects.append(dict(cohort=c,n=len(y),responders=int(y.sum()),nonresponders=int((1-y).sum()),full_AF35_g=a,residual_AF35_g=b,delta_g=b-a,retained_fraction=b/a,reference_full_g=expected,reference_abs_error=abs(a-expected)))
 seed=2026090801+k;rng=np.random.default_rng(seed);ir=np.flatnonzero(y==1);inn=np.flatnonzero(y==0)
 idx=np.concatenate([rng.choice(ir,size=(20000,len(ir)),replace=True),rng.choice(inn,size=(20000,len(inn)),replace=True)],axis=1)
 f.save(pd.DataFrame(idx,columns=[f'index_{i+1:03d}' for i in range(len(y))]).assign(draw=np.arange(1,20001)),f'inputs/{c}_PAIRED_BOOTSTRAP_INDICES_ZERO_BASED.tsv')
 ga=bhg(x,idx,len(ir));gb=bhg(res,idx,len(ir));delta=gb-ga;valid=np.isfinite(delta)
 boots.append(pd.DataFrame(dict(cohort=c,draw=np.arange(1,20001),full_g=ga,residual_g=gb,delta_g=delta,finite=valid,seed=seed)))
 q=np.quantile(delta[valid],[.025,.975]);summaries.append(dict(cohort=c,draws=20000,finite_draws=int(valid.sum()),nonfinite_draws=int((~valid).sum()),seed=seed,observed_delta_g=b-a,delta_q025=q[0],delta_q975=q[1],proportion_delta_lt_zero=float((delta[valid]<0).mean()),interpretation='DESCRIPTIVE_NOT_P_VALUE',residual_model_refit_in_bootstrap=False))
 cl=[]
 for omit,df in folds[folds.cohort==c].groupby('omitted_patient',sort=False):
  base=d.set_index('patient_id').loc[df.patient_id];yy=base.response.to_numpy();fx=base.AF35.to_numpy();rr=df.AF35_residual.to_numpy()
  aa=hg(fx,yy);bb=hg(rr,yy)
  cl.append(dict(cohort=c,omitted_patient=omit,n=len(yy),full_g=aa,residual_g=bb,delta_g=bb-aa,retained_fraction=bb/aa if aa!=0 else np.nan,refit_intercept=df.intercept.iloc[0],refit_broad_immune_beta=df.broad_immune_beta.iloc[0],residual_model_refit=True))
 loos+=cl;cd=pd.DataFrame(cl)
 looss.append(dict(cohort=c,folds=len(cd),positive_residual_g_n=int((cd.residual_g>0).sum()),residual_g_min=cd.residual_g.min(),residual_g_max=cd.residual_g.max(),retained_fraction_min=cd.retained_fraction.min(),retained_fraction_max=cd.retained_fraction.max()))
for s,x in [('CLINICAL_EFFECTS',effects),('PAIRED_BOOTSTRAP',pd.concat(boots)),('BOOTSTRAP_SUMMARY',summaries),('PATIENT_LOO',loos),('LOO_SUMMARY',looss)]:f.save(x,P+'_'+s+'.tsv')
f.save(pd.concat(joined),'inputs/PATIENT_ANALYSIS_VALUES_AFTER_FREEZE.tsv')
f.save(bindings,'QA/CLINICAL_SOURCE_BINDINGS.tsv')
npositive=sum(r['residual_AF35_g']>0 for r in effects)
pattern='RESIDUAL_POSITIVE_BOTH' if npositive==2 else ('RESIDUAL_POSITIVE_HETEROGENEOUS' if npositive==1 else ('RESIDUAL_NEAR_NULL_BOTH' if all(r['residual_AF35_g']==0 for r in effects) else 'RESIDUAL_REVERSED'))
f.js('QA/CLINICAL_COMPLETION_LOCK.json',dict(utc=f.now(),scientific_pattern=pattern,freeze_utc=lock['utc'],outcome_join_utc=(OUT/(P+'_OUTCOME_JOIN_TIMESTAMP.txt')).read_text().strip(),hashes={p.name:f.sha(p) for p in OUT.glob(P+'_*.tsv')}))
f.write(P+'_OUTCOME_FIREWALL_AUDIT.md','''# G 结局接入审计

来源审计先确定唯一历史broad_immune。PROTOCOL_PRECOMPUTATION_LOCK先于OLS；FREEZE_LOCK含原AF35、广义免疫值、全样本系数、残差和所有LOO重拟合分子值及哈希。03_join_and_analyze.py在再次核验所有输入与分子锁后才写OUTCOME_JOIN_TIMESTAMP并读取Analysis A患者response。标签文件哈希另列CLINICAL_SOURCE_BINDINGS。
临床分析没有再估计评分、改变基因或增减患者。响应只用于固定效应量比较和已预指定的分层配对bootstrap。LOO模型已在接入标签前拟合并锁定。分层bootstrap和仅使用冻结全样本残差的范围已在接入前写入协议。
本次为执行流程上的无结局残差构建，不声称历史项目或分析者从未见过总体效应。H主分析及其并列数值修正已完成才执行G接入；H不使用此处结果。
''')
print(pd.DataFrame(effects).to_string(index=False));print(pd.DataFrame(summaries).to_string(index=False));print(pd.DataFrame(looss).to_string(index=False));print(pattern)
