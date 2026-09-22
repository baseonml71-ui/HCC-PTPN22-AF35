from common import *
from scipy.stats import rankdata
lock=json.loads((OUT/(P+'_FREEZE_LOCK.json')).read_text(encoding='utf-8'))
assert not lock['outcome_labels_loaded'] and not (OUT/'QA/OUTCOME_JOIN_TIMESTAMP.txt').exists()
for p,h in lock['hashes'].items(): assert sha(OUT/p)==h,p
for p,h in lock['source_hashes'].items(): assert sha(p)==h,p
# Record the global lock hash before the first clinical source read.
join=now();write('QA/OUTCOME_JOIN_TIMESTAMP.txt',join)
js('QA/OUTCOME_OPEN_EVENT.json',dict(utc=join,freeze_utc=lock['utc'],freeze_lock_sha256=sha(OUT/(P+'_FREEZE_LOCK.json')),all_frozen_hashes_verified=True))
assert lock['utc']<join
sm=read(OUT/'AF35_GENE_NUMBER_SUBSET_MANIFEST.tsv');ss=read(OUT/'AF35_GENE_NUMBER_RESPONSE_BLIND_SCORES.tsv');ps=read(OUT/'AF35_STRUCTURED_REFERENCE_RESPONSE_BLIND_SCORES.tsv');stab=read(OUT/'AF35_GENE_NUMBER_SCORE_STABILITY.tsv')
sel=read(OUT/'AF35_STRUCTURED_REFERENCE_SELECTED_200.tsv').program_id.tolist()
expected={str(base/r.path):r.sha256 for base in [A,B] for r in read(base/'COMPLETE_SHA256_MANIFEST.tsv').itertuples()}
clinical_sources=[A/f'{c}_ANALYSIS_PATIENT_VALUES.tsv' for c in C]+[B/'AF35_MATCHED_NULL_PROGRAM_EFFECTS.tsv']
binding=[dict(path=str(p),sha256=sha(p),expected_sha256=expected[str(p)]) for p in clinical_sources]
assert all(r['sha256']==r['expected_sha256'] for r in binding);save(binding,'QA/CLINICAL_SOURCE_BINDINGS.tsv')
oldb=read(B/'AF35_MATCHED_NULL_PROGRAM_EFFECTS.tsv').set_index(['cohort','program_id'])
effects=[];jeffects=[];summ=[];labels=[];repro=[]
for c in C:
 ids=read(OUT/f'inputs/{c}_PATIENT_IDS.tsv').patient_id.tolist()
 lab=read(A/f'{c}_ANALYSIS_PATIENT_VALUES.tsv',usecols=['patient_id','response']).set_index('patient_id')
 assert lab.index.is_unique and set(lab.index)==set(ids) and set(lab.response)=={0,1}
 lab=lab.loc[ids];r=lab.response.to_numpy()==1
 labels.append(lab.reset_index().assign(cohort=c))
 x=ss[ss.cohort==c].pivot(index='patient_id',columns='subset_id',values='score').loc[ids,sm.subset_id].to_numpy()
 g=hg(x,r);af=g[-1]
 ed=pd.DataFrame(dict(cohort=c,subset_id=sm.subset_id,k=sm.k,n=len(ids),n_R=int(r.sum()),n_NR=int((~r).sum()),hedges_g=g));effects.append(ed)
 for k,d in ed.groupby('k'):
  v=d.hedges_g;sv=stab[(stab.cohort==c)&(stab.k==k)].spearman_rho
  summ.append(dict(cohort=c,k=k,subset_n=len(d),rho_median=sv.median(),rho_q025=sv.quantile(.025),rho_q975=sv.quantile(.975),rho_minimum=sv.min(),rho_maximum=sv.max(),g_median=v.median(),g_mean=v.mean(),g_q025=v.quantile(.025),g_q975=v.quantile(.975),g_positive_proportion=(v>0).mean(),g_minimum=v.min(),g_maximum=v.max(),full_AF35_reference_g=af))
 pdx=ps[ps.cohort==c].pivot(index='patient_id',columns='program_id',values='score').loc[ids]
 pg=hg(pdx.to_numpy(),r);pe=pd.DataFrame(dict(cohort=c,program_id=pdx.columns,n=len(ids),n_R=int(r.sum()),n_NR=int((~r).sum()),hedges_g=pg))
 err=np.max(np.abs(pg-oldb.loc[(c,pdx.columns),'hedges_g'].to_numpy()));assert err<1e-10
 repro.append(dict(cohort=c,program_n=len(pe),max_abs_error_vs_B=err,AF35_J1_J2_error=abs(af-pe.set_index('program_id').loc['AF35','hedges_g']),status='PASS'))
 jeffects.append(pe[pe.program_id.isin(['AF35']+sel)].copy())
save(pd.concat(labels),'inputs/RESPONSE_LABELS_AFTER_FREEZE.tsv');save(pd.concat(effects),'AF35_GENE_NUMBER_CLINICAL_EFFECTS.tsv');save(summ,'AF35_GENE_NUMBER_SUMMARY.tsv');save(pd.concat(jeffects),'AF35_STRUCTURED_REFERENCE_CLINICAL_EFFECTS.tsv');save(repro,'QA/CLINICAL_REPRODUCTION.tsv')
plane=pd.concat(jeffects).pivot(index='program_id',columns='cohort',values='hedges_g').loc[sel+['AF35']];save(plane.reset_index(),'AF35_STRUCTURED_REFERENCE_EFFECT_PLANE.tsv')
rp=plane.loc[sel];af=plane.loc['AF35'];ex=(rp>=af);counts={c:int(ex[c].sum()) for c in C};joint=int(ex.all(axis=1).sum())
cross=dict(selected_n=200,ERP_exceedance=counts[C[0]],GSE_exceedance=counts[C[1]],joint_exceedance=joint,ERP_AF35_g=af[C[0]],GSE_AF35_g=af[C[1]],effect_pearson=float(corr(rp[C[0]].to_numpy(),rp[C[1]].to_numpy())),effect_spearman=float(corr(rankdata(rp[C[0]]),rankdata(rp[C[1]]))),count_interpretation='DESCRIPTIVE_NOT_P_VALUE')
save([cross],'AF35_STRUCTURED_REFERENCE_CROSS_COHORT.tsv')
exceptional=[counts[c]<=lock['classification_max_exceedance'] for c in C]
state='COMPOSITION_SPECIFICITY_SUPPORTED' if all(exceptional) else ('MIXED_COHORT_RESULT' if any(exceptional) else 'AGGREGATION_AND_STRUCTURE_LARGELY_EXPLAIN_EFFECT')
interp=[]
for c in C:
 d=pd.DataFrame(summ).query('cohort == @c').sort_values('k');a=d.iloc[0];b=d[d.k==30].iloc[0]
 interp.append(dict(component='J1',cohort=c,state='DESCRIPTIVE_AGGREGATION_DIAGNOSTIC',rho_k5=a.rho_median,rho_k30=b.rho_median,rho_median_nondecreasing=bool((np.diff(d.rho_median)>=0).all()),g_k5=a.g_median,g_k30=b.g_median,g_median_nondecreasing=bool((np.diff(d.g_median)>=0).all()),g_range_width_k5=a.g_q975-a.g_q025,g_range_width_k30=b.g_q975-b.g_q025,boundary='聚合稳定性和效应分布诊断；不选择子集或最优基因数'))
interp.append(dict(component='J2',cohort='BOTH',state=state,boundary='固定描述性上5%规则；不是因果分解、正式置换检验或精确结构匹配；不能量化效应归因百分比'))
save(interp,P+'_INTERPRETATION.tsv')
write(P+'_OUTCOME_FIREWALL_AUDIT.md',f'''# J 结局防火墙

全局冻结：{lock['utc']}；首次临床表解析事件：{join}；二者严格先后。事件锁保存冻结JSON本身SHA256；全部{len(lock['hashes'])}个冻结文件及{len(lock['source_hashes'])}个来源哈希在接入前逐项校验通过。
01_freeze.py仅解析A/B/G/H明确无结局输入、方案/实现、来源清单与历史时间戳；两部分全部分数、六维结构量、距离、固定200身份均在最终全局锁内。02_clinical.py事件落盘后才解析A响应和B临床效应。项目保护仅二进制哈希，不等同读取结局表内容。检查历史评分代码时可见固定组人数；不声称分析者对全部项目历史结果从未知情。本轮不存在I输出读取或结果驱动参考选择。
结局接入后未改子集、评分、特征、K、选择或分类规则。第二语言R使用已锁定分子输入独立复算；技术修正须留痕，不允许改设计。
''')
js('QA/CLINICAL_COMPLETION_LOCK.json',dict(utc=now(),scientific_state=state,counts=counts,joint=joint,hashes={n:sha(OUT/n) for n in ['AF35_GENE_NUMBER_CLINICAL_EFFECTS.tsv','AF35_GENE_NUMBER_SUMMARY.tsv','AF35_STRUCTURED_REFERENCE_CLINICAL_EFFECTS.tsv','AF35_STRUCTURED_REFERENCE_CROSS_COHORT.tsv']}))
print(json.dumps(cross,indent=2));print(state)
