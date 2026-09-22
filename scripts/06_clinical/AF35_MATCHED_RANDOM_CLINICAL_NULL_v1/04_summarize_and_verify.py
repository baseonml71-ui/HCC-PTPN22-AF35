from pathlib import Path
import json,hashlib,datetime
import numpy as np
import pandas as pd

OUT=Path(__file__).resolve().parent
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write(d,n): d.to_csv(OUT/n,sep='\t',index=False)
def percentile(x,a): return 100*(np.sum(x<a)+.5*np.sum(x==a))/len(x)
lock=json.loads((OUT/'MATCHED_PROGRAM_FREEZE_LOCK.json').read_text())
for n,h in lock['hashes'].items(): assert sha(OUT/n)==h
join=datetime.datetime.fromisoformat((OUT/'OUTCOME_JOIN_TIMESTAMP.txt').read_text().strip().replace('Z','+00:00'))
assert datetime.datetime.fromisoformat(lock['time_utc'])<join
for row in pd.read_csv(OUT/'PATIENT_SCORE_FREEZE_SHA256.tsv',sep='\t').itertuples(): assert sha(OUT/row.path)==row.sha256
manifest=pd.read_csv(OUT/'AF35_MATCHED_NULL_PROGRAM_MANIFEST.tsv',sep='\t')
effects=pd.read_csv(OUT/'AF35_MATCHED_NULL_PROGRAM_EFFECTS.tsv',sep='\t')
patient=pd.read_csv(OUT/'AF35_MATCHED_NULL_PATIENT_PROGRAM_SCORES.tsv',sep='\t')
af=pd.read_csv(OUT/'AF35_ORDERED_GENES.tsv',sep='\t').gene.tolist()
sets={'AF35':af,**manifest.groupby('program_id',sort=True).matched_gene.apply(list).to_dict()}
verify=[]; summaries=[]
for cohort in ['ERP117672','GSE302495']:
    # Independent NumPy implementation starts from transformed gene expression, not R scores.
    x=pd.read_csv(OUT/f'{cohort}_SELECTED_TRANSFORMED_EXPRESSION.tsv',sep='\t',index_col=0)
    arr=x.to_numpy(); gz=(arr-arr.mean(axis=1,keepdims=True))/arr.std(axis=1,ddof=1,keepdims=True)
    ids=x.columns.to_list(); index={g:i for i,g in enumerate(x.index)}
    meta=pd.read_csv(OUT/'ANALYSIS_PATIENT_ELIGIBILITY.tsv',sep='\t').query('cohort==@cohort').set_index('patient_id').loc[ids]
    included=meta.analysis_included.to_numpy(); response=meta.response_group.eq('Responder').to_numpy()
    ep=effects.query('cohort==@cohort').set_index('program_id')
    sp=patient.query('cohort==@cohort').pivot(index='program_id',columns='patient_id',values='score')[ids]
    for pid,genes in sets.items():
        vals=gz[[index[g] for g in genes]].mean(axis=0); vals=(vals-vals.mean())/vals.std(ddof=1)
        r=vals[included & response]; nr=vals[included & ~response]; n=len(r)+len(nr)
        g=(1-3/(4*(n-2)-1))*(r.mean()-nr.mean())/np.sqrt(((len(r)-1)*r.var(ddof=1)+(len(nr)-1)*nr.var(ddof=1))/(n-2))
        se=np.max(abs(vals-sp.loc[pid].to_numpy())); ge=abs(g-ep.loc[pid,'hedges_g'])
        verify.append(dict(cohort=cohort,program_id=pid,max_patient_score_absolute_error=se,hedges_g_absolute_error=ge,pass_check=bool(se<1e-10 and ge<1e-10)))
    random=ep.drop('AF35').hedges_g.to_numpy(); a=float(ep.loc['AF35','hedges_g'])
    q=np.quantile(random,[.025,.25,.5,.75,.975],method='linear')
    summaries.append(dict(cohort=cohort,n=int(ep.loc['AF35','n']),n_R=int(ep.loc['AF35','n_R']),n_NR=int(ep.loc['AF35','n_NR']),program_n=2000,AF35_g=a,random_q025=q[0],random_q25=q[1],random_median_g=q[2],random_q75=q[3],random_q975=q[4],AF35_percentile=percentile(random,a),n_random_ge_AF35=int(np.sum(random>=a)),fraction_random_ge_AF35=np.mean(random>=a),n_random_positive=int(np.sum(random>0)),fraction_random_positive=np.mean(random>0)))
verification=pd.DataFrame(verify); write(verification,'INDEPENDENT_NUMERIC_VERIFICATION.tsv'); assert verification.pass_check.all()
summary=pd.DataFrame(summaries); write(summary,'AF35_MATCHED_NULL_PRIMARY_SUMMARY.tsv')
cross=effects.pivot(index='program_id',columns='cohort',values='hedges_g').rename(columns={'ERP117672':'g_ERP','GSE302495':'g_GSE'})
cross['mean_g']=cross[['g_ERP','g_GSE']].mean(axis=1); cross['min_g']=cross[['g_ERP','g_GSE']].min(axis=1)
cross['both_positive']=cross.g_ERP.gt(0)&cross.g_GSE.gt(0)
cross['both_ge_corresponding_AF35']=cross.g_ERP.ge(cross.loc['AF35','g_ERP'])&cross.g_GSE.ge(cross.loc['AF35','g_GSE'])
write(cross.reset_index(),'AF35_MATCHED_NULL_CROSS_COHORT.tsv')
random=cross.drop('AF35'); a=cross.loc['AF35']; joint=[]
for metric in ['mean_g','min_g']:
    q=np.quantile(random[metric],[.025,.25,.5,.75,.975],method='linear')
    joint.append(dict(metric=metric,AF35_value=a[metric],AF35_percentile=percentile(random[metric].to_numpy(),a[metric]),random_q025=q[0],random_q25=q[1],random_median=q[2],random_q75=q[3],random_q975=q[4],n_random_ge_AF35=int((random[metric]>=a[metric]).sum()),fraction_random_ge_AF35=(random[metric]>=a[metric]).mean(),n_both_positive=int(random.both_positive.sum()),fraction_both_positive=random.both_positive.mean(),n_both_ge_AF35=int(random.both_ge_corresponding_AF35.sum()),fraction_both_ge_AF35=random.both_ge_corresponding_AF35.mean(),program_n=2000))
joint=pd.DataFrame(joint); write(joint,'AF35_MATCHED_NULL_JOINT_SUMMARY.tsv')
(OUT/'AF35_MATCHED_NULL_OUTCOME_FIREWALL_AUDIT.md').write_text(f'''# 结局防火墙和程序冻结审计

AF35_MATCHED_NULL_OUTCOME_FIREWALL_GATE = PASS
AF35_MATCHED_NULL_MANIFEST_FREEZE_GATE = PASS

共同程序冻结时间 UTC：{lock['time_utc']}。第一次响应文件解析前的时间标记：{join.isoformat()}。
程序 SHA256：{lock['hashes']['AF35_MATCHED_NULL_PROGRAM_MANIFEST.tsv']}。
生成程序仅打开基因测量特征、基因定义和无结局评分权威锁；不读取患者标签。六维匹配参数和全部程序先行冻结。评分脚本先核验全部锁，计算并哈希 78×2001 份患者分数，再解析两份原始响应映射；ERP 的 5 例不可评估患者只参与冻结分数标准化，分析时排除。
冻结哈希和来源锁重验通过。完整 AF35 的 g 复现通过；独立 NumPy 从 transformed expression 重算 4002 个程序-队列结果，分数和效应误差均 <1e-10。
这是代码与文件级执行顺序审计，不是操作系统级访问隔离证明。分析者已从附件知道历史 AF35 总体效应；不能称分析者完全盲法，也不能据此恢复原始 AF35 选择防火墙。
''',encoding='utf-8')
print(summary.to_string(index=False));print(joint.to_string(index=False));print('INDEPENDENT_NUMERIC_CHECKS=4002/4002 PASS')
