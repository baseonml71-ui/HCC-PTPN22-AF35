from common import *
import shutil, platform
from scipy.stats import rankdata
for d in ['inputs','QA','figures','authority']: (OUT/d).mkdir(exist_ok=True)
assert not (OUT/(P+'_FREEZE_LOCK.json')).exists()
pass  # Non-computational private authorization copy excluded.
if not (OUT/'QA/PROTECTED_BEFORE.tsv').exists():
 paths=protected_paths();print('Hashing protected files',len(paths),sum(p.stat().st_size for p in paths)/1e9,'GB',flush=True)
 save([dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in paths],'QA/PROTECTED_BEFORE.tsv')
 js('QA/INITIALIZATION_LOCK.json',dict(utc=now(),analysis_I_concurrently_active_per_author=True,excluded_active_directory=str(ROOT/'AF35_FUNCTIONAL_ROLE_STATE_ARCHITECTURE_v1'),excluded_dependencies=['node_modules','.pnpm','__pycache__','.git'],scope='All other project files including external authoritative inputs; hash only, no clinical parsing'))
expected={}
for base in [A,B,G,H]:
 for r in read(base/'COMPLETE_SHA256_MANIFEST.tsv').itertuples():expected[str(base/r.path)]=r.sha256
sources=[A/'AF35_ORDERED_GENES.tsv',B/'AF35_ORDERED_GENES.tsv',B/'AF35_MATCHED_NULL_PROGRAM_MANIFEST.tsv',B/'MATCHED_PROGRAM_FREEZE_LOCK.json',B/'PATIENT_SCORE_FREEZE_SHA256.tsv',B/'OUTCOME_JOIN_TIMESTAMP.txt',B/'AF35_MATCHED_NULL_PATIENT_PROGRAM_SCORES_RESPONSE_BLIND.tsv',G/'AF35_BROAD_IMMUNE_FREEZE_LOCK.json',G/'AF35_BROAD_IMMUNE_RESPONSE_BLIND_PATIENT_VALUES.tsv',G/'inputs/BROAD_IMMUNE_GENES.tsv',G/'scripts/02_freeze_residuals.py',H/'scripts/01_freeze_and_analyze.py',H/'AF35_PROGRAM_COHERENCE_RANDOM_REFERENCE_RESULTS.tsv',H/'AF35_PROGRAM_COHERENCE_WITHIN_COHORT_SUMMARY.tsv']
for c in C:sources.extend([B/f'{c}_SELECTED_TRANSFORMED_EXPRESSION.tsv',A/f'{c}_GENE_Z_AND_FULL_SCORES_RESPONSE_BLIND.tsv',H/f'inputs/{c}_PATIENT_IDS.tsv',H/f'inputs/{c}_GENE_Z_RESPONSE_BLIND.tsv',G/f'inputs/{c}_PATIENT_IDS.tsv',G/f'inputs/{c}_AF35_BROAD_TRANSFORMED_EXPRESSION.tsv'])
rows=[dict(path=str(p),sha256=sha(p),expected_sha256=expected[str(p)],status='PASS' if sha(p)==expected[str(p)] else 'FAIL') for p in sources]
save(rows,'AF35_AGGREGATION_SOURCE_BINDING_AUDIT.tsv');assert all(r['status']=='PASS' for r in rows)
genes=read(A/'AF35_ORDERED_GENES.tsv').gene.tolist()
assert genes==read(B/'AF35_ORDERED_GENES.tsv').gene.tolist() and len(set(genes))==35 and 'PTPN22' not in genes
shutil.copyfile(A/'AF35_ORDERED_GENES.tsv',OUT/'inputs/AF35_ORDERED_GENES.tsv')
manifest=read(B/'AF35_MATCHED_NULL_PROGRAM_MANIFEST.tsv');bl=json.loads((B/'MATCHED_PROGRAM_FREEZE_LOCK.json').read_text())
bt=(B/'OUTCOME_JOIN_TIMESTAMP.txt').read_text().strip()
assert not bl['outcome_labels_loaded'] and bl['time_utc']<bt
assert sha(B/'AF35_MATCHED_NULL_PROGRAM_MANIFEST.tsv')==bl['hashes']['AF35_MATCHED_NULL_PROGRAM_MANIFEST.tsv']
for rr in read(B/'PATIENT_SCORE_FREEZE_SHA256.tsv').itertuples():assert sha(B/rr.path)==rr.sha256
gl=json.loads((G/'AF35_BROAD_IMMUNE_FREEZE_LOCK.json').read_text());assert not gl['outcome_labels_loaded']
assert sha(G/'AF35_BROAD_IMMUNE_RESPONSE_BLIND_PATIENT_VALUES.tsv')==gl['hashes']['AF35_BROAD_IMMUNE_RESPONSE_BLIND_PATIENT_VALUES.tsv']
sets={pid:d.sort_values('position').matched_gene.tolist() for pid,d in manifest.groupby('program_id',sort=True)}
assert len(sets)==2000 and len({tuple(sorted(gs)) for gs in sets.values()})==2000
assert all(len(gs)==len(set(gs))==35 for gs in sets.values())
assert all(d.sort_values('position').AF35_gene.tolist()==genes for _,d in manifest.groupby('program_id'))
shutil.copyfile(B/'AF35_MATCHED_NULL_PROGRAM_MANIFEST.tsv',OUT/'inputs/EXACT_ANALYSIS_B_RANDOM_PROGRAM_MANIFEST.tsv')
protocol='''# Analysis J：执行前方案冻结

仅本独立模块写入。A/B/G/H为直接权威；不读取、等待或写入Analysis I。用户声明I并行活动，其完整输出目录从不可变审计排除；其余项目文件（依赖缓存除外）仅二进制哈希，保护外部权威输入，不解析其中结局内容。
固定无符号等权AF35原35基因及顺序，排除PTPN22；PTPN22仍为BIOLOGICAL_ANCHOR。ERP35、GSE38患者与G/H完全一致，保留原ERP40/GSE38的基因z标准化空间。每个程序和子集先平均gene-z，再按原40/38人作评分z，与A/B完全相同；同时保存未二次z的等权均值。取35/38人后不重标准化评分。该正仿射变换不改变rho或g。
J1：k=5,10,15,20,25,30每组2000次独立均匀无放回选基因，NumPy PCG64 seed=2026090809；每子集内基因唯一，子集间允许重复（独立抽样规则，不因重复重抽）。原序排序；两队列共用一个清单。k35仅固定全AF35一次。无结局下计算全部患者分数及Spearman/描述性Pearson；average ranks。每k报告中位数、type7/linear分位数2.5/97.5%、极值。临床阶段加g均值、中位数、分位范围、正值比例和全AF35参考。
J2：B原始2000个35基因程序逐字节复用，原ID/成员不变。仅六维：每队列median pairwise Spearman、PC1 variance explained、程序评分与G固定BroadImmune的Spearman。相关矩阵为average ranks后Pearson，上三角595对中位数；H实现PCA在实际35/38人内重新中心化、按样本SD标准化、SVD平方奇异值比例，不构造PC1临床评分。BroadImmune固定11基因PTPRC,HLA-DRA,LST1,TYROBP,FCER1G,CD74,MS4A1,CD79A,CD3D,CD3E,NKG7；原40/38人gene-z均值，不另作评分z。先与G原值复现再使用G精确冻结值。
各维按全部2000随机程序均值/样本SD(ddof=1)转z；同变换用于AF35；欧氏距离；距离升序、精确距离并列时program_id字典序，固定最近200，禁止其他K。匹配QC标准差差值=(AF35-selected median)/all2000 SD，保留符号；报告选中距离median/max及未选minimum。最近参照不等于精确匹配或生物等价。
严格顺序：分子输入→J1清单→两部分全部无结局评分→六维结构量→距离与固定200身份→全部哈希锁→结局接入。临床表/响应标签在最终锁前不得读取。历史脚本可检查公式；其中固定组人数属于方法元数据，不能宣称分析者历史上从未见过项目结局背景。
临床只接A冻结ANALYSIS_PATIENT_VALUES.tsv的patient_id/response；g沿用A/B，R-NR，合并组内样本SD，J=1-3/[4(n-2)-1]。与B原效应逐一核对。所有子集/程序只作扰动诊断，患者为生物学重复，不由子集计数计算P值，不运行ROC/AUC/分类器。
描述分类预先操作化：单队列“明显上尾”为selected中g>=AF35者<=10/200，即描述性上5%。两队列均满足（joint自动<=10）为COMPOSITION_SPECIFICITY_SUPPORTED；仅一队列满足为MIXED_COHORT_RESULT；两队列均>10为AGGREGATION_AND_STRUCTURE_LARGELY_EXPLAIN_EFFECT。最后标签只表示不满足稀少超越规则，不能据此量化“解释大部分”或证明因果分解。连续效应、具体计数和残余结构差异优先于标签；不可改变阈值。联合超越为同一程序同时满足两队列>=条件，原数值直接比较；不叫P值。
J1不拟合阈值或强制单调。聚合稳定性、效应范围缩窄、效应中位数变化分别解释；不把35称最优，不比较挑选最佳子集，不宣称各基因必需、唯一预测、ICI特异性或因果。
独立R从B变换表达重新计算gene-z、全部J1评分/相关/g/汇总、全部J2六维/标准化/距离/200身份/效应/超越/跨队列相关；数值容差1e-10，身份和计数严格一致。Python独立重放随机种子。技术PASS与科学方向分开。图180mm可编辑SVG、单页PDF、600dpi PNG，渲染审核。结束STOP_FOR_AUTHOR_SCIENTIFIC_REVIEW，不改手稿或Figures1–6，不开始K。
'''
write(P+'_PROTOCOL_FREEZE_v1.md',protocol)
js('QA/PROTOCOL_PRECOMPUTATION_LOCK.json',dict(utc=now(),outcome_labels_loaded=False,sha256=sha(OUT/(P+'_PROTOCOL_FREEZE_v1.md')),seed=2026090809,classification_max_exceedance=10,K=200))
rng=np.random.default_rng(2026090809);subs=[]
for k in [5,10,15,20,25,30]:
 for j in range(2000):
  q=sorted(rng.choice(35,k,replace=False).tolist())
  subs.append(dict(subset_id=f'K{k:02d}_{j+1:04d}',k=k,draw=j+1,positions=';'.join(str(i+1) for i in q),genes=';'.join(genes[i] for i in q),seed=2026090809))
subs.append(dict(subset_id='AF35',k=35,draw=1,positions=';'.join(str(i+1) for i in range(35)),genes=';'.join(genes),seed=2026090809))
sm=pd.DataFrame(subs);save(sm,'AF35_GENE_NUMBER_SUBSET_MANIFEST.tsv')
w=np.zeros((35,len(sm)))
for j,row in sm.iterrows():w[[int(i)-1 for i in row.positions.split(';')],j]=1/row.k
broad=read(G/'AF35_BROAD_IMMUNE_RESPONSE_BLIND_PATIENT_VALUES.tsv')
bgenes=read(G/'inputs/BROAD_IMMUNE_GENES.tsv').gene.tolist()
assert len(bgenes)==11 and set(bgenes)==set('PTPRC HLA-DRA LST1 TYROBP FCER1G CD74 MS4A1 CD79A CD3D CD3E NKG7'.split())
shutil.copyfile(G/'inputs/BROAD_IMMUNE_GENES.tsv',OUT/'inputs/BROAD_IMMUNE_GENES.tsv')
features=[];stability=[];j1scores=[];j2scores=[];recover=[]
features_by_id={pid:dict(program_id=pid) for pid in ['AF35']+list(sets)}
for c in C:
 ids=read(H/f'inputs/{c}_PATIENT_IDS.tsv').patient_id.tolist()
 assert ids==read(G/f'inputs/{c}_PATIENT_IDS.tsv').patient_id.tolist() and len(set(ids))==len(ids)==(35 if c==C[0] else 38)
 save({'patient_id':ids},f'inputs/{c}_PATIENT_IDS.tsv')
 raw=read(B/f'{c}_SELECTED_TRANSFORMED_EXPRESSION.tsv',index_col=0)
 z=raw.sub(raw.mean(axis=1),axis=0).div(raw.std(axis=1,ddof=1),axis=0)
 x=z[ids].T;olda=read(A/f'{c}_GENE_Z_AND_FULL_SCORES_RESPONSE_BLIND.tsv').set_index('patient_id')
 assert np.max(np.abs(x[genes].to_numpy()-olda.loc[ids,genes].to_numpy()))<1e-12
 hx=read(H/f'inputs/{c}_GENE_Z_RESPONSE_BLIND.tsv').set_index('patient_id')
 assert x.columns.tolist()==hx.columns.tolist() and np.max(np.abs(x.to_numpy()-hx.loc[ids].to_numpy()))<1e-12
 save(x.rename_axis('patient_id').reset_index(),f'inputs/{c}_GENE_Z_RESPONSE_BLIND.tsv')
 brraw=read(G/f'inputs/{c}_AF35_BROAD_TRANSFORMED_EXPRESSION.tsv',index_col=0)
 brz=brraw.sub(brraw.mean(axis=1),axis=0).div(brraw.std(axis=1,ddof=1),axis=0)
 br=broad[broad.cohort==c].set_index('patient_id').loc[ids,'BroadImmune'].to_numpy()
 err=np.max(np.abs(brz.loc[bgenes,ids].mean(axis=0).to_numpy()-br));assert err<1e-12
 save({'patient_id':ids,'BroadImmune':br},f'inputs/{c}_BROAD_IMMUNE_RESPONSE_BLIND.tsv')
 sr=z.loc[genes].T.to_numpy()@w;sz=(sr-sr.mean(axis=0))/sr.std(axis=0,ddof=1);ix=[raw.columns.get_loc(i) for i in ids]
 sr=sr[ix];sz=sz[ix];full=olda.loc[ids,'AF35'].to_numpy()
 assert np.max(np.abs(sz[:,-1]-full))<1e-12
 rho=corr(rankdata(sz,axis=0),rankdata(full)[:,None]);pear=corr(sz,full[:,None])
 stability.append(pd.DataFrame({'cohort':c,'subset_id':sm.subset_id,'k':sm.k,'spearman_rho':rho,'pearson_r':pear}))
 j1scores.append(pd.DataFrame({'cohort':c,'patient_id':np.repeat(ids,len(sm)),'subset_id':np.tile(sm.subset_id,len(ids)),'k':np.tile(sm.k,len(ids)),'equal_weight_mean_gene_z':sr.ravel(),'score':sz.ravel()}))
 for pid,gs in [('AF35',genes)]+list(sets.items()):
  xx=x[gs].to_numpy();rr=rankdata(xx,axis=0);cm=np.corrcoef(rr,rowvar=False);vv=cm[np.triu_indices(35,1)]
  zz=(xx-xx.mean(axis=0))/xx.std(axis=0,ddof=1);ss=np.linalg.svd(zz,compute_uv=False)
  pr=z.loc[gs].mean(axis=0);ps=(pr-pr.mean())/pr.std(ddof=1);pv=ps.loc[ids].to_numpy()
  f=features_by_id[pid];f[c+'__median_pairwise']=float(np.median(vv));f[c+'__PC1_variance']=float(ss[0]**2/(ss*ss).sum());f[c+'__immune_spearman']=float(corr(rankdata(pv),rankdata(br)))
  j2scores.append(pd.DataFrame({'cohort':c,'patient_id':ids,'program_id':pid,'score':pv}))
 recover.append(dict(cohort=c,n=len(ids),standardization_n=len(raw.columns),broad_max_error=err,AF35_score_max_error=float(np.max(np.abs(sz[:,-1]-full))),status='PASS'))
 print(c,'all response-blind scores and features complete',flush=True)
save(pd.concat(j1scores),'AF35_GENE_NUMBER_RESPONSE_BLIND_SCORES.tsv');save(pd.concat(stability),'AF35_GENE_NUMBER_SCORE_STABILITY.tsv');save(pd.concat(j2scores),'AF35_STRUCTURED_REFERENCE_RESPONSE_BLIND_SCORES.tsv');save(recover,'QA/MOLECULAR_RECOVERY.tsv')
ft=pd.DataFrame(features_by_id.values()).set_index('program_id');fn=ft.columns.tolist();rd=ft.loc[list(sets)];mu=rd.mean();sd=rd.std(ddof=1);assert (sd>0).all()
zt=(ft-mu)/sd;dist=np.sqrt(((zt.loc[list(sets)]-zt.loc['AF35'])**2).sum(axis=1))
ds=pd.DataFrame({'program_id':dist.index,'STRUCTURAL_DISTANCE_TO_AF35':dist.values}).sort_values(['STRUCTURAL_DISTANCE_TO_AF35','program_id']).reset_index(drop=True);ds['distance_rank']=np.arange(1,2001);ds['selected']=ds.distance_rank<=200
sel=ds.head(200);selected=sel.program_id.tolist()
save(ft.reset_index(),'AF35_STRUCTURED_REFERENCE_FEATURES.tsv');save(zt.reset_index(),'AF35_STRUCTURED_REFERENCE_STANDARDIZED_FEATURES.tsv');save(pd.DataFrame({'feature':fn,'random_program_mean':mu.values,'random_program_SD':sd.values}),'AF35_STRUCTURED_REFERENCE_STANDARDIZATION.tsv');save(ds,'AF35_STRUCTURAL_DISTANCE_TO_AF35.tsv');save(sel,'AF35_STRUCTURED_REFERENCE_SELECTED_200.tsv')
qc=[]
for col in fn:
 a=ft.loc[selected,col];q=a.quantile([.025,.5,.975]);qc.append(dict(feature=col,AF35_value=ft.loc['AF35',col],selected_median=q.loc[.5],selected_q025=q.loc[.025],selected_q975=q.loc[.975],all2000_median=rd[col].median(),all2000_q025=rd[col].quantile(.025),all2000_q975=rd[col].quantile(.975),standardized_difference=(ft.loc['AF35',col]-q.loc[.5])/sd[col]))
save(qc,'AF35_STRUCTURED_REFERENCE_MATCHING_QC.tsv');save([dict(K=200,median_selected_distance=sel.STRUCTURAL_DISTANCE_TO_AF35.median(),maximum_selected_distance=sel.STRUCTURAL_DISTANCE_TO_AF35.max(),minimum_nonselected_distance=ds.iloc[200].STRUCTURAL_DISTANCE_TO_AF35)],'AF35_STRUCTURED_REFERENCE_DISTANCE_QC.tsv')
oldh=read(H/'AF35_PROGRAM_COHERENCE_RANDOM_REFERENCE_RESULTS.tsv').set_index('program_id')
for c in C:
 for key in ['median_pairwise','PC1_variance']:assert np.max(np.abs(ft.loc[list(sets),c+'__'+key]-oldh.loc[list(sets),c+'__'+key]))<1e-10
write('AF35_STRUCTURED_REFERENCE_PROVENANCE_AUDIT.md',f'''# J2 来源核验

B原始清单逐字节复制，SHA256={sha(B/'AF35_MATCHED_NULL_PROGRAM_MANIFEST.tsv')}；2000个精确ID、35个唯一成员、原序AF35位置映射均通过。B冻结{bl['time_utc']}，结局接入{bt}；原锁outcome_labels_loaded=false。表达及原评分哈希与B患者分数锁一致。G免疫值与G结局前冻结哈希一致，11基因公式重新复现。H相关与PCA实现原样复算；2000程序两队列结果与H数值一致。未使用I注释，未生成新参考程序或免疫基因宇宙。
''')
paths=[p for p in OUT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='PROTECTED_BEFORE.tsv']
lock=dict(utc=now(),outcome_labels_loaded=False,seed=2026090809,K=200,subset_n=12000,full_program_n=1,random_program_n=2000,classification_max_exceedance=10,python=platform.python_version(),numpy=np.__version__,hashes={str(p.relative_to(OUT)):sha(p) for p in sorted(paths)},source_hashes={r['path']:r['sha256'] for r in rows})
js(P+'_FREEZE_LOCK.json',lock)
print('J GLOBAL FREEZE COMPLETE',lock['utc'],'No clinical tables or response labels parsed',flush=True)
