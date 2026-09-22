import importlib, json, platform
import numpy as np
import pandas as pd
from scipy.stats import rankdata
f=importlib.import_module('00_initialize');OUT=f.OUT;P=f.PREFIX;B=f.B
def put(x,s):f.save(x,P+'_'+s+'.tsv')
assert not (OUT/(P+'_FREEZE_LOCK.json')).exists(),'Existing freeze; do not overwrite'
genes=f.read(OUT/'inputs/AF35_ORDERED_GENES.tsv').gene.tolist()
manifest=f.read(B/'AF35_MATCHED_NULL_PROGRAM_MANIFEST.tsv')
bl=json.loads((B/'MATCHED_PROGRAM_FREEZE_LOCK.json').read_text())
joined=(B/'OUTCOME_JOIN_TIMESTAMP.txt').read_text().strip()
assert bl['outcome_labels_loaded'] is False and f.sha(B/'AF35_MATCHED_NULL_PROGRAM_MANIFEST.tsv')==bl['hashes']['AF35_MATCHED_NULL_PROGRAM_MANIFEST.tsv']
assert bl['time_utc']<joined
sets={k:d.sort_values('position').matched_gene.tolist() for k,d in manifest.groupby('program_id',sort=True)}
assert len(sets)==2000 and len({tuple(sorted(v)) for v in sets.values()})==2000
assert all(len(v)==len(set(v))==35 for v in sets.values())
assert all(d.sort_values('position').AF35_gene.tolist()==genes for _,d in manifest.groupby('program_id'))
f.save(manifest,'inputs/EXACT_ANALYSIS_B_RANDOM_PROGRAM_MANIFEST.tsv')
f.write(P+'_RANDOM_REFERENCE_PROVENANCE.md',f'''# 随机程序来源

只复用 Analysis B 已冻结的 2,000 个可测性匹配、各 35 基因程序，逐行/逐程序完全相同，未生成新随机程序。
来源：{B/'AF35_MATCHED_NULL_PROGRAM_MANIFEST.tsv'}。
SHA256：{f.sha(B/'AF35_MATCHED_NULL_PROGRAM_MANIFEST.tsv')}。
原锁时间 {bl['time_utc']}；原结局接入时间 {joined}。前者早于后者，锁内 outcome_labels_loaded=false。
本模块未读取随机程序临床效应、临床超越结果或响应标签；仅读取清单、来源哈希和无结局表达矩阵。
该分布为可测性匹配随机程序参照，非正式置换零分布；原匹配来自 40/38 人评分空间，本模块结构分析使用冻结的 35/38 人患者集合，匹配程序本身不变。
''')
seed=2026090808;rng=np.random.default_rng(seed);parts=[];seen=set();attempts=0
while len(parts)<5000:
 q=tuple(sorted(rng.choice(35,17,replace=False).tolist()));attempts+=1
 if q not in seen:seen.add(q);parts.append(q)
spl=pd.DataFrame([dict(partition=i+1,half_A_positions=';'.join(str(v+1) for v in q),half_A_genes=';'.join(genes[v] for v in q),half_B_genes=';'.join(g for j,g in enumerate(genes) if j not in q),seed=seed) for i,q in enumerate(parts)])
f.save(spl,'inputs/SPLIT_HALF_PARTITION_MANIFEST.tsv')
protocol='''# H：响应盲法结构分析执行前冻结

仅本目录交付，PTPN22 仍为生物锚点；AF35 为原序、无符号、等权 35 基因状态读出，不替代评分。患者为重复。
主集合：从 Analysis E 无结局分子输入读取 patient_id，ERP117672=35，GSE302495=38；顺序固定。不读取响应或效应表。
表达来自 Analysis B 已冻结 SELECTED_TRANSFORMED_EXPRESSION.tsv：ERP log2(TPM+1)、GSE log2(CPM+1)。先在原评分空间 ERP40/GSE38 按样本 SD 标准化每个基因，再按固定 ID 取 ERP35/GSE38，须与 Analysis A 每基因 z 值逐项一致。此历史基因 z 用于等权 gene-to-rest 与 split-half。PCA 单独在实际 35/38 人中再次中心化并以样本 SD 标准化，只保存解释方差，不生成临床 PC1 分数。
Spearman=average ranks 后 Pearson，35×35 完整矩阵；595 个唯一无对角配对按原序 i=1..34、j=i+1..35 输出。报告 mean/median、正值比例、type7 2.5/97.5%分位数。每基因与其余34基因均值 Spearman；连接度为与其余34基因 rho 的均值。
Split-half：NumPy PCG64 seed=2026090808，5000个唯一17/18分割，两队列共用；17基因半边作为规范标识，排除重复，互补18基因不再次算一种分割。报告median、type7 2.5/97.5%范围和2r/(1+r)校正的中位数；r=-1时校正无定义并显式记录，不选择最优分割。
跨队列：595配对向量 Pearson/Spearman、35连接度向量 Pearson/Spearman，均描述性；不把595对当独立生物重复，不计算其显著性。
原样复用 Analysis B 2000程序，无重抽程序；每程序各队列 median pairwise rho、PC1 variance、median gene-to-rest；跨队列两种 covariance 和两种 connectivity 相关。参照 median、2.5/97.5%范围、>=AF35计数及严格小于AF35的比例；不是P值。随机程序每个再跑5000分割需2000万半分计算，本次省略该可选项，保留全部必需参照指标。
GSE104580为可选，本次不纳入，也不搜索或重建其输入；主推断仅上述两队列。
预定义描述标签：两队列 median pairwise>0、median gene-to-rest>0且median split-half>0为内部正向一致；跨队列两种covariance相关>0为方向保留。如果内部一致且两种cross相关与两队列median pairwise都超过各自随机参照97.5%分位数，COHERENT_AND_CROSS_COHORT_CONSERVED；内部一致、cross均>0但未满足全部参照条件，COHERENT_BUT_NOT_EXCEPTIONAL；内部一致而cross任一<=0，WITHIN_COHORT_COHERENT_CROSS_COHORT_WEAK；否则NOT_STRUCTURALLY_COHERENT。这些只是预先声明的操作标签，无显著性或生物阈值含义；以连续量和参照位置为主要解释。
第二语言R从冻结表达矩阵独立复算所有核心指标、2000随机程序和5000分割；数值容差1e-10，身份/顺序严格一致。数值/技术PASS不以正向结果为条件。独立180mm图、PDF渲染、600dpi PNG及可编辑文字SVG。完成即作者科学审阅停止，不改旧分析、手稿或主图。
'''
f.write(P+'_PROTOCOL_FREEZE_v1.md',protocol)
lockfiles=[OUT/(P+'_PROTOCOL_FREEZE_v1.md'),OUT/'inputs/SPLIT_HALF_PARTITION_MANIFEST.tsv',OUT/'inputs/EXACT_ANALYSIS_B_RANDOM_PROGRAM_MANIFEST.tsv',OUT/'inputs/AF35_ORDERED_GENES.tsv',OUT/'scripts/01_freeze_and_analyze.py',OUT/(P+'_SOURCE_BINDING_AUDIT.tsv')]
lockfiles += [OUT/f'inputs/{c}_PATIENT_IDS.tsv' for c in f.COHORTS]
f.js(P+'_FREEZE_LOCK.json',dict(utc=f.now(),outcome_labels_loaded=False,seed=seed,partitions=5000,partition_attempts=attempts,source_hashes={r.path:r.sha256 for r in f.read(OUT/(P+'_SOURCE_BINDING_AUDIT.tsv')).itertuples()},hashes={str(p.relative_to(OUT)):f.sha(p) for p in lockfiles}))
print('H PROTOCOL FROZEN; starting response-blind metrics',flush=True)
def corrcols(x,y):
 x=x-x.mean(axis=0);y=y-y.mean(axis=0)
 return (x*y).sum(axis=0)/np.sqrt((x*x).sum(axis=0)*(y*y).sum(axis=0))
def paircorr(x,y):return float(corrcols(x,y))
def cross(a,b):
 return dict(covariance_pearson=paircorr(a['vector'],b['vector']),covariance_spearman=paircorr(rankdata(a['vector']),rankdata(b['vector'])),connectivity_pearson=paircorr(a['connectivity'],b['connectivity']),connectivity_spearman=paircorr(rankdata(a['connectivity']),rankdata(b['connectivity'])))
ii,jj=np.triu_indices(35,1)
def metrics(x):
 ranks=rankdata(x,axis=0);m=np.corrcoef(ranks,rowvar=False);v=m[ii,jj]
 z=(x-x.mean(axis=0))/x.std(axis=0,ddof=1);s=np.linalg.svd(z,compute_uv=False);ev=s*s/(s*s).sum()
 rest=(x.sum(axis=1,keepdims=True)-x)/34
 gr=corrcols(ranks,rankdata(rest,axis=0))
 return dict(matrix=m,vector=v,connectivity=(m.sum(axis=1)-1)/34,rest=gr,median_pairwise=float(np.median(v)),mean_pairwise=float(v.mean()),positive_pair_fraction=float((v>0).mean()),pair_q025=float(np.quantile(v,.025)),pair_q975=float(np.quantile(v,.975)),PC1_variance=float(ev[0]),PC2_variance=float(ev[1]),median_gene_to_rest=float(np.median(gr)),min_gene_to_rest=float(gr.min()),max_gene_to_rest=float(gr.max()))
data={};index={};af={};within=[];pc=[];grs=[];splits=[]
w=np.zeros((35,5000))
for k,q in enumerate(parts):w[list(q),k]=1
for c in f.COHORTS:
 raw=f.read(B/f'{c}_SELECTED_TRANSFORMED_EXPRESSION.tsv',index_col=0)
 assert raw.index.is_unique and np.isfinite(raw.to_numpy()).all()
 z=(raw.sub(raw.mean(axis=1),axis=0)).div(raw.std(axis=1,ddof=1),axis=0)
 ids=f.read(OUT/f'inputs/{c}_PATIENT_IDS.tsv').patient_id.tolist()
 x=z[ids].T;assert set(genes).issubset(x.columns)
 old=f.read(f.A/f'{c}_GENE_Z_AND_FULL_SCORES_RESPONSE_BLIND.tsv').set_index('patient_id')
 assert np.max(np.abs(x[genes].to_numpy()-old.loc[ids,genes].to_numpy()))<1e-12
 f.save(x.rename_axis('patient_id').reset_index(),'inputs/'+c+'_GENE_Z_RESPONSE_BLIND.tsv')
 data[c]=x.to_numpy();index[c]={g:i for i,g in enumerate(x.columns)}
 xa=x[genes].to_numpy();m=metrics(xa);af[c]=m
 put(pd.DataFrame(m['matrix'],columns=genes).assign(gene=genes)[['gene']+genes],'GENE_CORRELATION_'+c)
 within.append(dict(cohort=c,n=len(ids),gene_n=35,pair_n=595,**{k:v for k,v in m.items() if np.isscalar(v)}))
 pc.extend([dict(cohort=c,component=i+1,variance_explained=m['PC1_variance' if i==0 else 'PC2_variance'],diagnostic_only=True) for i in range(2)])
 grs.extend([dict(cohort=c,gene_rank=i+1,gene=g,rho=m['rest'][i],connectivity=m['connectivity'][i]) for i,g in enumerate(genes)])
 ra=xa@w/17;rb=xa@(1-w)/18;rho=corrcols(rankdata(ra,axis=0),rankdata(rb,axis=0))
 sb=np.divide(2*rho,1+rho,out=np.full_like(rho,np.nan),where=rho!=-1)
 splits.extend([dict(cohort=c,partition=i+1,rho=r,spearman_brown=sb[i]) for i,r in enumerate(rho)])
 within[-1].update(split_half_median=float(np.median(rho)),split_half_q025=float(np.quantile(rho,.025)),split_half_q975=float(np.quantile(rho,.975)),spearman_brown_median=float(np.nanmedian(sb)))
 print(c,'AF35 input and metrics complete',flush=True)
put(within,'WITHIN_COHORT_SUMMARY');put(pc,'PCA');put(grs,'GENE_TO_REST');put(splits,'SPLIT_HALF')
ac=cross(*[af[c] for c in f.COHORTS]);put([dict(cohort_pair='ERP117672__GSE302495',pair_n=595,gene_n=35,**ac)],'CROSS_COHORT')
put([dict(gene_i=genes[i],gene_j=genes[j],ERP117672=af[f.COHORTS[0]]['matrix'][i,j],GSE302495=af[f.COHORTS[1]]['matrix'][i,j]) for i,j in zip(ii,jj)],'COVARIANCE_VECTORS')
random=[]
for k,(program,gs) in enumerate(sets.items()):
 mm={c:metrics(data[c][:,[index[c][g] for g in gs]]) for c in f.COHORTS}
 row=dict(program_id=program)
 for c in f.COHORTS:
  for metric in ['median_pairwise','PC1_variance','median_gene_to_rest']:row[c+'__'+metric]=mm[c][metric]
 row.update(cross(*[mm[c] for c in f.COHORTS]));random.append(row)
 if (k+1)%500==0:print('H matched programs',k+1,flush=True)
rd=pd.DataFrame(random);put(rd,'RANDOM_REFERENCE_RESULTS')
obs={c+'__'+m:af[c][m] for c in f.COHORTS for m in ['median_pairwise','PC1_variance','median_gene_to_rest']};obs.update(ac)
rr=[]
for metric,o in obs.items():
 v=rd[metric].to_numpy();q=np.quantile(v,[.025,.5,.975])
 rr.append(dict(metric=metric,AF35_observed=o,reference_n=2000,reference_median=q[1],reference_q025=q[0],reference_q975=q[2],n_ge_AF35=int((v>=o).sum()),empirical_position_fraction_strictly_below=float((v<o).mean()),interpretation='DESCRIPTIVE_REFERENCE_NOT_P_VALUE'))
put(rr,'RANDOM_REFERENCE_SUMMARY')
positive=all(r['median_pairwise']>0 and r['median_gene_to_rest']>0 and r['split_half_median']>0 for r in within)
cr=ac['covariance_pearson']>0 and ac['covariance_spearman']>0
keys=['covariance_pearson','covariance_spearman']+[c+'__median_pairwise' for c in f.COHORTS]
exceptional=all(r['AF35_observed']>r['reference_q975'] for r in rr if r['metric'] in keys)
pattern='NOT_STRUCTURALLY_COHERENT' if not positive else ('WITHIN_COHORT_COHERENT_CROSS_COHORT_WEAK' if not cr else ('COHERENT_AND_CROSS_COHORT_CONSERVED' if exceptional else 'COHERENT_BUT_NOT_EXCEPTIONAL'))
f.js('QA/PRIMARY_COMPLETION_LOCK.json',dict(utc=f.now(),scientific_pattern=pattern,outcome_labels_loaded=False,random_programs=2000,split_half_partitions=5000,python=platform.python_version(),numpy=np.__version__,hashes={p.name:f.sha(p) for p in OUT.glob(P+'_*.tsv')}))
f.write(P+'_OUTCOME_FIREWALL_AUDIT.md','''# H 结局防火墙审计

本轮 H 主分析进程无响应标签或临床效应表读入。此前项目背景及用户附件包含历史总体效应，不声称分析者历史上从未见过效应；本次执行隔离以输入与顺序为证。
脚本只允许无结局表达、原序基因、原随机清单、既有无结局patient_id及哈希/时间锁。读取B的 OUTCOME_JOIN_TIMESTAMP.txt 只检查原清单冻结先后；未读B临床结果。保护文件仅二进制哈希，无结局内容解析。
00_initialize.py 建立独立保护基线；01_freeze_and_analyze.py 先冻结协议、分割清单、来源哈希，再计算H，PRIMARY_COMPLETION_LOCK.json固定结果。G标签接入在H主分析结束以后执行。所有图和解释来自结构量，无响应、效应量、临床超越或结果驱动选择。
GSE104580及随机程序split-half可选项未运行；未增加新队列或程序。见协议中的预先选择。
''')
print('H PRIMARY COMPLETE',pattern,flush=True)
