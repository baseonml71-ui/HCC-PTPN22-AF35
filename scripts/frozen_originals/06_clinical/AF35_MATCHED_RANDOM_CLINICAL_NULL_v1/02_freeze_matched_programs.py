from pathlib import Path
import hashlib,json,datetime,sys,re
import numpy as np
import pandas as pd

OUT=Path(__file__).resolve().parent
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write(df,name): df.to_csv(OUT/name,sep='\t',index=False)
assert not (OUT/'MATCHED_PROGRAM_FREEZE_LOCK.json').exists(), 'Already frozen'
for row in pd.read_csv(OUT/'SCORING_AUTHORITY_LOCK.tsv',sep='\t').itertuples(): assert sha(OUT/row.path)==row.sha256
d=pd.read_csv(OUT/'SOURCE_GENE_FEATURES_RESPONSE_BLIND.tsv',sep='\t')
genes=pd.read_csv(OUT/'AF35_ORDERED_GENES.tsv',sep='\t').gene.tolist()
cohorts=['ERP117672','GSE302495']
wide=d.pivot(index='gene',columns='cohort',values=['mean_expression','sample_sd','detection_fraction','variable'])
wide.columns=[f'{b}_{a}' for a,b in wide.columns]
valid=np.ones(len(wide),bool)
for c in cohorts: valid &= wide[f'{c}_variable'].fillna(False).astype(bool).to_numpy()
technical=wide.index.str.match(r'^(MT-|RPL|RPS)')
excluded=wide.index.isin(genes+['PTPN22'])
eligible=valid & ~technical & ~excluded
wide['is_AF35']=wide.index.isin(genes)
wide['eligible_candidate']=eligible
wide['exclusion_reason']=np.select([~valid,excluded,technical],['NOT_FINITE_VARIABLE_BOTH','AF35_OR_PTPN22','PROJECT_TECHNICAL_PREFIX'],default='ELIGIBLE')
features=[f'{c}_{f}' for c in cohorts for f in ['mean_expression','sample_sd','detection_fraction']]
candidate=wide.loc[eligible,features].astype(float).sort_index()
target=wide.loc[genes,features].astype(float)
assert len(candidate)>=200 and np.isfinite(candidate).all().all() and np.isfinite(target).all().all()
center=candidate.mean(); scale=candidate.std(ddof=1); assert (scale>0).all()
cx=((candidate-center)/scale).to_numpy(); tx=((target-center)/scale).to_numpy()
for f in features: wide[f+'_matching_z']=(wide[f].astype(float)-center[f])/scale[f]
write(wide.reset_index(),'AF35_MATCHED_NULL_MATCHING_FEATURES.tsv')
write(pd.DataFrame(dict(feature=features,center=center.values,sample_sd=scale.values,n_candidates=len(candidate))),'MATCHING_STANDARDIZATION_PARAMETERS.tsv')
names=candidate.index.to_numpy(); pools=[]; distances=[]; poolrows=[]
for i,g in enumerate(genes):
    dist=np.sqrt(((cx-tx[i])**2).sum(axis=1))
    order=np.argsort(dist,kind='stable')[:200]
    pools.append(order); distances.append(dist)
    poolrows.extend(dict(position=i+1,AF35_gene=g,pool_rank=r+1,matched_gene=names[k],distance=dist[k]) for r,k in enumerate(order))
write(pd.DataFrame(poolrows),'MATCHED_NEIGHBOR_POOLS.tsv')
rng=np.random.default_rng(20260908); rows=[]; selected=np.empty((2000,35),dtype=int)
for j in range(2000):
    used=set()
    for step,i in enumerate(rng.permutation(35)):
        available=[int(k) for k in pools[i] if int(k) not in used]
        assert len(available)>=166
        k=int(rng.choice(available)); used.add(k); selected[j,i]=k
        rows.append(dict(program_id=f'RANDOM_{j+1:04d}',position=int(i+1),AF35_gene=genes[i],matched_gene=names[k],distance=distances[i][k],pool_rank=int(np.flatnonzero(pools[i]==k)[0]+1),sampling_step=step+1,master_seed=20260908))
manifest=pd.DataFrame(rows).sort_values(['program_id','position'])
assert len(manifest)==70000 and manifest.groupby('program_id').matched_gene.nunique().eq(35).all()
assert not manifest.matched_gene.isin(genes+['PTPN22']).any()
write(manifest,'AF35_MATCHED_NULL_PROGRAM_MANIFEST.tsv')
delta=cx[selected].mean(axis=1)-tx.mean(axis=0)
qc=pd.DataFrame(dict(feature=features,AF35_mean_z=tx.mean(axis=0),random_mean_z=cx[selected].mean(axis=(0,1)),mean_bias=delta.mean(axis=0),q95_absolute_program_mean_bias=np.quantile(abs(delta),.95,axis=0)))
qc['pass']=qc.mean_bias.abs().le(.25) & qc.q95_absolute_program_mean_bias.le(.5)
write(qc,'MATCHING_BALANCE_QC.tsv')
pos=[]
for i,g in enumerate(genes):
    dd=np.array([distances[i][k] for k in selected[:,i]])
    pos.append(dict(position=i+1,AF35_gene=g,nearest_distance=distances[i][pools[i][0]],pool_max_distance=distances[i][pools[i][-1]],sampled_median_distance=np.median(dd),sampled_q95_distance=np.quantile(dd,.95),unique_matched_genes=len(set(selected[:,i])),**{f'{f}_mean_bias':float((cx[selected[:,i]]-tx[i]).mean(axis=0)[k]) for k,f in enumerate(features)}))
write(pd.DataFrame(pos),'POSITION_MATCHING_QC.tsv')
unique_sets=len({tuple(sorted(names[x])) for x in selected})
status='PASS' if qc['pass'].all() else 'FAIL'
(OUT/'AF35_MATCHED_NULL_MATCHING_QC.md').write_text(f'''# 匹配可行性核验

AF35_MATCHED_NULL_MATCHING_FEASIBILITY_GATE = {status}

两个矩阵均有限且变异的共同基因：{int(valid.sum())}；排除后共同候选：{len(candidate)}。完整逐基因纳入/排除原因在 MATCHING_FEATURES。
六维来源和固定参数见方案；两个队列检测零点有效且有变化，纳入自然检测比例。AF35 多数基因检测率为 1，存在饱和。
35 个位置均有 200 最近候选；生成 2000 个程序、70000 位置，程序内恰好 35 独特基因。独特无序程序集合数：{unique_sets}。实际涉及 {len(np.unique(selected))} 个候选基因。
逐维平均偏差 <=0.25 且 95% 程序绝对均值偏差 <=0.50：{status}。详见 MATCHING_BALANCE_QC.tsv；这仅保证预设全程序测量特征平衡，不保证逐基因完美匹配或共表达/功能相同。每位置最差匹配在 POSITION_MATCHING_QC.tsv 全部保留。
本阶段未读取响应标签或随机程序效应，不因结局调整池、种子、特征或程序。
''',encoding='utf-8')
assert status=='PASS','MATCHING_FEASIBILITY_FAIL: STOP before outcomes'
frozen=['AF35_MATCHED_NULL_PROTOCOL_FREEZE_v1.md','02_freeze_matched_programs.py','AF35_ORDERED_GENES.tsv','AF35_MATCHED_NULL_PROGRAM_MANIFEST.tsv','AF35_MATCHED_NULL_MATCHING_FEATURES.tsv','MATCHING_STANDARDIZATION_PARAMETERS.tsv','MATCHED_NEIGHBOR_POOLS.tsv','MATCHING_BALANCE_QC.tsv','POSITION_MATCHING_QC.tsv','AF35_MATCHED_NULL_MATCHING_QC.md','SCORING_AUTHORITY_LOCK.tsv']
lock=dict(time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),master_seed=20260908,rng='NumPy default_rng PCG64',numpy_version=np.__version__,python=sys.executable,program_n=2000,genes_per_program=35,candidate_n=len(candidate),unique_gene_n=len(np.unique(selected)),unique_program_sets=unique_sets,outcome_labels_loaded=False,hashes={n:sha(OUT/n) for n in frozen})
(OUT/'MATCHED_PROGRAM_FREEZE_LOCK.json').write_text(json.dumps(lock,indent=2),encoding='utf-8')
write(pd.DataFrame([dict(path=n,sha256=h) for n,h in lock['hashes'].items()]),'MATCHED_PROGRAM_MANIFEST_SHA256.tsv')
(OUT/'AF35_MATCHED_NULL_PROGRAM_MANIFEST.sha256').write_text(sha(OUT/'AF35_MATCHED_NULL_PROGRAM_MANIFEST.tsv')+'  AF35_MATCHED_NULL_PROGRAM_MANIFEST.tsv\n')
print(qc.to_string(index=False)); print(json.dumps({k:v for k,v in lock.items() if k!='hashes'},indent=2))
