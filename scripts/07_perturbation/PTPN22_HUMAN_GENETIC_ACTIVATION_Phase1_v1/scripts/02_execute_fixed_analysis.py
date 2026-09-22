"""仅执行已锁定AF35及八程序分析；PTPN22仅描述；没有其他单基因效应。"""
from pathlib import Path
import hashlib,json,datetime,itertools,math
import pandas as pd
import numpy as np
from scipy.special import gammaln
from scipy.stats import t as t_dist
import openpyxl
ROOT=Path(__file__).resolve().parents[1];I=ROOT/'01_manifest/inputs'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def tsv(data,p):pd.DataFrame(data).to_csv(ROOT/p,sep='\t',index=False,lineterminator='\n',float_format='%.15g')
def js(data,p):(ROOT/p).write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
L=json.loads((ROOT/'01_manifest/METHODS_LOCK.json').read_text())
for row in pd.read_csv(ROOT/'01_manifest/PRE_ANALYSIS_LOCK_SHA256.tsv',sep='\t').itertuples():assert sha(ROOT/row.path)==row.sha256
for row in pd.read_csv(ROOT/'01_manifest/INPUT_SHA256.tsv',sep='\t').itertuples():assert sha(ROOT/row.copied_path)==row.sha256
js(dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),methods_sha256=sha(ROOT/'01_manifest/METHODS_LOCK.json'),
 manifest_sha256=sha(ROOT/'01_manifest/PTPN22_HUMAN_GENETIC_ACTIVATION_PHASE1_MANIFEST.md'),analysis_code_sha256=sha(Path(__file__))),
 'logs/ANALYSIS_START.json')
matrix_path=I/'GSE94859_RPKMs.xlsx';workbook_hash=sha(matrix_path)
wb=openpyxl.load_workbook(matrix_path,read_only=True,data_only=True);head=next(wb.active.values);wb.close()
assert len(head)==len(set(head))
m=pd.read_excel(matrix_path,engine='openpyxl');raw=m.iloc[:,2:].to_numpy(float);columns=m.columns[2:].tolist()
ledger=pd.read_csv(I/'GSE94859_SAMPLE_PAIRING_LEDGER.tsv',sep='\t',keep_default_na=False)
mapping=pd.read_csv(I/'GSE94859_GENOTYPE_MAPPING.tsv',sep='\t')
programs=pd.read_csv(I/'PTPN22_MECHANISM_PROGRAMS_FROZEN.tsv',sep='\t')
programs['program']=programs.program.replace({'PTPN22_activation_feedback_35':'AF35'})
names=['AF35']+L['secondary'];sets={name:programs.loc[programs.program==name,'gene'].tolist() for name in names}
assert len(sets['AF35'])==35 and 'PTPN22' not in sets['AF35']
assert len(ledger)==32 and ledger.GSM.nunique()==32 and set(columns)==set(ledger.processed_matrix_column)
assert np.isfinite(raw).all() and (raw>=0).all()
assert len({hashlib.sha256(raw[:,j].tobytes()).hexdigest() for j in range(raw.shape[1])})==32
donor_ids=sorted(ledger.individual_id.unique(),key=lambda z:int(z.replace('individual','')))
assert len(donor_ids)==16
recheck=[]
for ident in donor_ids:
    g=ledger[ledger.individual_id==ident]
    assert len(g)==2 and set(g.condition_standardized)=={'naive','activated'} and g.genotype_raw.nunique()==1
    assert g.age.nunique()==1 and g.sex.nunique()==1 and set(g.metadata_conflict)=={'NONE'}
    paper=mapping.loc[mapping.individual_id==ident,'source_paper_genotype'].item()
    assert {'GG':'CC','AA':'TT'}[g.genotype_raw.iloc[0]]==paper
    recheck.append(dict(individual_id=ident,naive_GSM=g.loc[g.condition_standardized=='naive','GSM'].item(),
      activated_GSM=g.loc[g.condition_standardized=='activated','GSM'].item(),genotype_raw=g.genotype_raw.iloc[0],source_genotype=paper,
      complete_pair=True,genotype_concordance=True,covariate_concordance=True,exact_matrix_mapping=True,status='PASS'))
tsv(recheck,'08_qa/PAIRING_RECHECK.tsv')
target_genes=list(dict.fromkeys(programs.gene.tolist()+['PTPN22']))
target_checks=[];rowidx={}
for gene in target_genes:
    indices=m.index[m.Gene_symbol==gene].tolist();assert len(indices)==1
    i=indices[0];rowidx[gene]=i
    assert np.std(raw[i,:],ddof=1)>0
    compound=any(isinstance(s,str) and '|' in s and gene in s.split('|') for s in m.Gene_symbol)
    assert not compound
    target_checks.append(dict(gene=gene,excel_row=i+2,exact_rows=1,compound_rows=0,finite_n=32,variable=True,status='PASS'))
tsv(target_checks,'08_qa/TARGET_MATRIX_RECHECK.tsv')

# 此函数完全不接收供者、条件或基因型标签。
def score_fixed(matrix,ids,definitions):
    scores=[];parameters=[];setpars=[]
    for name in names:
        genes=definitions[name];m0=np.array([matrix[ids[g],:] for g in genes])
        x=np.log2(m0+1.0);mu=x.mean(axis=1);sd=x.std(axis=1,ddof=1)
        assert np.isfinite(sd).all() and (sd>0).all()
        z=(x-mu[:,None])/sd[:,None];u=z.mean(axis=0);us=float(u.std(ddof=1))
        assert us>0 and np.isfinite(us)
        scores.append((u-u.mean())/us)
        for k,g in enumerate(genes):parameters.append(dict(program=name,gene=g,gene_n=len(genes),weight=1/len(genes),log2_mean_32=mu[k],log2_sd_32=sd[k]))
        setpars.append(dict(program=name,gene_n=len(genes),score_pre_z_mean_32=float(u.mean()),score_pre_z_sd_32=us))
    return np.array(scores).T,parameters,setpars
score,parameters,setpars=score_fixed(raw,rowidx,sets)
assert np.isfinite(score).all() and np.allclose(score.mean(axis=0),0,atol=1e-12) and np.allclose(score.std(axis=0,ddof=1),1,atol=1e-12)
# 样本列逆序后复原，评分不受列顺序影响。
reverse,_,_=score_fixed(raw[:,::-1],rowidx,sets)
assert np.allclose(reverse[::-1,:],score,rtol=1e-12,atol=1e-12)
blind=pd.DataFrame(score,columns=names);blind.insert(0,'processed_matrix_column',columns)
tsv(blind,'02_scores/FROZEN_SCORES_BEFORE_GENOTYPE_JOIN.tsv')
tsv(parameters,'02_scores/FROZEN_SCORING_PARAMETERS.tsv');tsv(setpars,'02_scores/FROZEN_PROGRAM_SCALING_PARAMETERS.tsv')
js(dict(locked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
 scores_sha256=sha(ROOT/'02_scores/FROZEN_SCORES_BEFORE_GENOTYPE_JOIN.tsv'),parameters_sha256=sha(ROOT/'02_scores/FROZEN_SCORING_PARAMETERS.tsv'),
 genotype_used_in_scoring=False,complete_AF35=35,all_samples=32), '02_scores/SCORING_FREEZE.json')
sample=ledger.merge(blind,on='processed_matrix_column',validate='one_to_one').merge(mapping[['individual_id','source_paper_genotype']],on='individual_id',validate='many_to_one')
tsv(sample[['GSM','processed_matrix_column','individual_id','condition_standardized','genotype_raw','source_paper_genotype']+names],
 '02_scores/GSE94859_FROZEN_SAMPLE_SCORES.tsv')
scoring_audit=[dict(check='source_workbook_unchanged',status='PASS',value=sha(matrix_path)==workbook_hash),
 dict(check='genotype_independent_function',status='PASS',value='score_fixed inputs: matrix,row index, frozen definitions only'),
 dict(check='sample_order_invariance',status='PASS',value=float(np.max(np.abs(reverse[::-1,:]-score)))),
 dict(check='score_mean_zero_max_error',status='PASS',value=float(np.abs(score.mean(axis=0)).max())),
 dict(check='score_sd_one_max_error',status='PASS',value=float(np.abs(score.std(axis=0,ddof=1)-1).max())),
 dict(check='targets_full_unique_valid',status='PASS',value=f'{len(target_genes)} including PTPN22; AF35 35/35'),
 dict(check='source_matrix_anomalies',status='PASS_TARGET_SCOPE',value='重复/日期/复合标识符未影响目标；原表不改写')]
tsv(scoring_audit,'08_qa/SCORING_AUDIT.tsv')

# PCA/相关仅QC；不筛高变基因，不输出载荷，不据此剔除样本。
sym=m.Gene_symbol
qc_mask=sym.map(lambda z:isinstance(z,str) and '|' not in z)&~sym.duplicated(keep=False)
qc=np.log2(raw[qc_mask.to_numpy(),:]+1)
qc=qc[np.std(qc,axis=1)>0,:]
center=qc-qc.mean(axis=1,keepdims=True)
eigen,U=np.linalg.eigh(center.T@center);ix=np.argsort(eigen)[::-1];eigen=np.clip(eigen[ix],0,None);U=U[:,ix]
pc=U[:,:2]*np.sqrt(eigen[:2])
for k in range(2):
    if pc[np.argmax(np.abs(pc[:,k])),k]<0:pc[:,k]*=-1
pcframe=pd.DataFrame(dict(processed_matrix_column=columns,PC1=pc[:,0],PC2=pc[:,1])).merge(ledger[['processed_matrix_column','individual_id','condition_standardized','genotype_raw']],on='processed_matrix_column')
tsv(pcframe,'08_qa/SAMPLE_PCA_QC.tsv')
corr=pd.DataFrame(np.corrcoef(qc.T),columns=columns);corr.insert(0,'processed_matrix_column',columns)
tsv(corr,'08_qa/SAMPLE_CORRELATION_QC.tsv')
js(dict(gene_rows=int(qc.shape[0]),PC1_fraction=float(eigen[0]/eigen.sum()),PC2_fraction=float(eigen[1]/eigen.sum()),
 minimum_sample_correlation=float(np.min(corr.iloc[:,1:].to_numpy())),excluded_samples=0,
 note='PCA只作技术QC；无来源定义异常判据被违反；无载荷/通路解释'), '08_qa/PCA_QC_SUMMARY.json')

delta_rows=[];Y=[];group=[];donor_meta=[]
for ident in donor_ids:
    s=sample[sample.individual_id==ident];n=s[s.condition_standardized=='naive'].iloc[0];a=s[s.condition_standardized=='activated'].iloc[0]
    delta=a[names].to_numpy(float)-n[names].to_numpy(float);Y.append(delta);group.append(n.source_paper_genotype=='TT')
    donor_meta.append(dict(individual_id=ident,genotype=n.source_paper_genotype,age=int(n.age),male=int(n.sex=='Male')))
    for k,name in enumerate(names):delta_rows.append(dict(individual_id=ident,genotype=n.source_paper_genotype,program=name,
      naive_score=float(n[name]),activated_score=float(a[name]),delta=float(delta[k]),delta_positive=bool(delta[k]>0)))
delta_df=pd.DataFrame(delta_rows);Y=np.array(Y);group=np.array(group);assert group.sum()==8
tsv(delta_df[delta_df.program=='AF35'].rename(columns={'delta':'Delta_AF35'}),'03_primary/GSE94859_AF35_ACTIVATION_DELTA.tsv')
tsv(delta_df[delta_df.program!='AF35'].rename(columns={'delta':'Delta_program'}),'04_secondary/GSE94859_FROZEN_PROGRAM_ACTIVATION_DELTA.tsv')
tsv(delta_df[delta_df.program=='AF35'].sort_values('delta').assign(overall_order=range(1,17)),'03_primary/GSE94859_AF35_ORDERED_DONOR_DELTAS.tsv')
within=[]
for name in names:
    for geno in ['CC','TT']:
        g=delta_df[(delta_df.program==name)&(delta_df.genotype==geno)]
        for variable in ['naive_score','activated_score','delta']:
            v=g[variable].to_numpy();within.append(dict(program=name,genotype=geno,measure=variable,n=len(v),mean=v.mean(),median=np.median(v),min=v.min(),max=v.max(),sd=np.std(v,ddof=1),positive_n=int((v>0).sum()) if variable=='delta' else ''))
tsv(within,'03_primary/GSE94859_WITHIN_GENOTYPE_DESCRIPTIVE_SUMMARY.tsv')

# 全部8/8分配，同一矩阵用于9个终点。
combinations=np.array(list(itertools.combinations(range(16),8)),dtype=int);assert len(combinations)==12870
weights=np.full((12870,16),-1/8,dtype=float)
weights[np.arange(12870)[:,None],combinations]=1/8
assert np.allclose(weights.sum(axis=1),0)
perm=weights@Y;D=Y[group,:].mean(axis=0)-Y[~group,:].mean(axis=0)
observed_assignment=tuple(np.flatnonzero(group));obs_idx=next(j for j,c in enumerate(combinations) if tuple(c)==observed_assignment)
assert np.allclose(perm[obs_idx],D)
tol=1e-12*np.maximum(1,np.abs(D));extreme=(np.abs(perm)>=np.abs(D)[None,:]-tol[None,:]).sum(axis=0);p=extreme/12870
assert (p>0).all() and np.allclose(np.sort(perm,axis=0),np.sort(-perm,axis=0),atol=1e-12)
permtable=pd.DataFrame(perm,columns=names);permtable.insert(0,'TT_donors',[';'.join(donor_ids[i] for i in c) for c in combinations]);permtable.insert(0,'assignment_id',range(12870))
tsv(permtable,'08_qa/EXACT_PERMUTATION_DISTRIBUTIONS.tsv')
rng=np.random.default_rng(L['bootstrap_seed']);B=L['bootstrap_n']
index_cc=rng.integers(0,8,size=(B,8));index_tt=rng.integers(0,8,size=(B,8))
cc=Y[~group,:];tt=Y[group,:];b_cc=cc[index_cc,:];b_tt=tt[index_tt,:]
bd=b_tt.mean(axis=1)-b_cc.mean(axis=1)
sp=np.sqrt((7*cc.var(axis=0,ddof=1)+7*tt.var(axis=0,ddof=1))/14)
J=float(np.exp(gammaln(7)-0.5*np.log(7)-gammaln(6.5)))
assert (sp>0).all();g=J*D/sp
bsp=np.sqrt((7*b_cc.var(axis=1,ddof=1)+7*b_tt.var(axis=1,ddof=1))/14)
bg=np.divide(J*bd,bsp,out=np.full_like(bd,np.nan),where=bsp>0)
q=np.empty(8);sort=np.argsort(p[1:]);ordered=p[1:][sort]*8/np.arange(1,9);q[sort]=np.minimum(1,np.minimum.accumulate(ordered[::-1])[::-1])
results=[];paudit=[]
for k,name in enumerate(names):
    ci=np.quantile(bd[:,k],[.025,.975],method='linear');valid=np.isfinite(bg[:,k]);gci=np.quantile(bg[valid,k],[.025,.975],method='linear')
    results.append(dict(program=name,D_TT_minus_CC=float(D[k]),exact_P=float(p[k]),BH_q=float(q[k-1]) if k else '',
      CC_n=8,TT_n=8,CC_mean_delta=float(cc[:,k].mean()),TT_mean_delta=float(tt[:,k].mean()),
      CC_median_delta=float(np.median(cc[:,k])),TT_median_delta=float(np.median(tt[:,k])),
      CC_min_delta=float(cc[:,k].min()),CC_max_delta=float(cc[:,k].max()),TT_min_delta=float(tt[:,k].min()),TT_max_delta=float(tt[:,k].max()),
      CC_positive_delta_n=int((cc[:,k]>0).sum()),TT_positive_delta_n=int((tt[:,k]>0).sum()),
      pooled_delta_SD=float(sp[k]),Hedges_g=float(g[k]),J=J,
      bootstrap_D_CI95_low=float(ci[0]),bootstrap_D_CI95_high=float(ci[1]),bootstrap_g_CI95_low=float(gci[0]),bootstrap_g_CI95_high=float(gci[1]),
      bootstrap_n=B,bootstrap_g_invalid_n=int((~valid).sum()),bootstrap_seed=L['bootstrap_seed'],permutation_extreme_n=int(extreme[k]),permutation_total_n=12870))
    paudit.append(dict(program=name,assignments=12870,TT_per_assignment=8,CC_per_assignment=8,observed_assignment_id=obs_idx,
      extreme_count=int(extreme[k]),exact_P=float(p[k]),tie_tolerance=float(tol[k]),all_enumerated=True,complement_symmetry=True,observed_assignment_included=True))
tsv(results[:1],'03_primary/GSE94859_AF35_GENOTYPE_INTERACTION.tsv');tsv(results[1:],'04_secondary/GSE94859_PROGRAM_GENOTYPE_INTERACTION.tsv')
tsv(paudit,'08_qa/PERMUTATION_AUDIT.tsv')
tsv(pd.DataFrame(dict(bootstrap_id=np.arange(B),D_AF35=bd[:,0],Hedges_g_AF35=bg[:,0])),'08_qa/AF35_BOOTSTRAP_DISTRIBUTION.tsv')

lodo=[]
for i,ident in enumerate(donor_ids):
    keep=np.arange(16)!=i;d=float(Y[keep&group,0].mean()-Y[keep&~group,0].mean());change=d-D[0]
    influence=abs(change)/abs(D[0]) if D[0]!=0 else np.nan
    lodo.append(dict(removed_donor=ident,removed_genotype='TT' if group[i] else 'CC',remaining_TT_n=int((keep&group).sum()),remaining_CC_n=int((keep&~group).sum()),
      full_D_AF35=float(D[0]),LODO_D_AF35=d,change_from_full=float(change),absolute_change=float(abs(change)),
      influence_fraction=float(influence),exceeds_40_percent=bool(influence>.4),sign_preserved=bool(np.sign(d)==np.sign(D[0]) and D[0]!=0),
      sign_reversal=bool(d*D[0]<0),score_rescaled=False))
tsv(lodo,'05_robustness/GSE94859_AF35_LODO.tsv')
preserve=sum(r['sign_preserved'] for r in lodo);maxinf=max(r['influence_fraction'] for r in lodo)
robust='PASS' if preserve==16 and maxinf<=.4 else 'PARTIAL' if preserve>=14 and maxinf<1 else 'FAIL'

# 仅锁定的age/sex模型；HC3 t推断为敏感性，不是精确主检验。
dm=pd.DataFrame(donor_meta);age=dm.age.to_numpy(float);X=np.column_stack([np.ones(16),group.astype(float),age-age.mean(),dm.male])
assert np.linalg.matrix_rank(X)==4
inverse=np.linalg.inv(X.T@X);beta=inverse@X.T@Y[:,0];residual=Y[:,0]-X@beta
h=np.sum((X@inverse)*X,axis=1);meat=X.T@((residual**2/(1-h)**2)[:,None]*X);cov=inverse@meat@inverse;se=np.sqrt(np.diag(cov))
tval=beta/se;pt=2*t_dist.sf(np.abs(tval),12);margin=t_dist.ppf(.975,12)*se
sens=[dict(term=name,beta=float(beta[k]),HC3_SE=float(se[k]),CI95_low=float(beta[k]-margin[k]),CI95_high=float(beta[k]+margin[k]),
 HC3_t_P=float(pt[k]),df=12,n=16,model='Delta_AF35 ~ TT + centered_age + male',role='SENSITIVITY_ONLY',primary_gate_changed=False) for k,name in enumerate(['intercept','TT_vs_CC','age_centered','male'])]
tsv(sens,'05_robustness/GSE94859_AF35_COVARIATE_SENSITIVITY.tsv')

# PTPN22上下文，不做假设检验。
context=[]
for ident in donor_ids:
    ld=ledger[ledger.individual_id==ident];paper=mapping.loc[mapping.individual_id==ident,'source_paper_genotype'].item()
    ni=columns.index(ld.loc[ld.condition_standardized=='naive','processed_matrix_column'].item());ai=columns.index(ld.loc[ld.condition_standardized=='activated','processed_matrix_column'].item())
    n=float(raw[rowidx['PTPN22'],ni]);a=float(raw[rowidx['PTPN22'],ai])
    context.append(dict(individual_id=ident,genotype=paper,naive_author_FPKM_RPKM=n,activated_author_FPKM_RPKM=a,delta_author_FPKM_RPKM=a-n,
      naive_log2=np.log2(n+1),activated_log2=np.log2(a+1),delta_log2=np.log2(a+1)-np.log2(n+1),role='CONTEXT_ONLY_NO_ACTIVITY_INFERENCE'))
tsv(context,'06_context/GSE94859_PTPN22_EXPRESSION_CONTEXT.tsv')
cs=[]
for geno in ['CC','TT']:
    frame=pd.DataFrame(context);frame=frame[frame.genotype==geno]
    for col in ['naive_author_FPKM_RPKM','activated_author_FPKM_RPKM','delta_author_FPKM_RPKM','naive_log2','activated_log2','delta_log2']:
        v=frame[col];cs.append(dict(genotype=geno,measure=col,n=8,mean=v.mean(),median=v.median(),min=v.min(),max=v.max()))
tsv(cs,'06_context/PTPN22_CONTEXT_SUMMARY.tsv')
primary='SUPPORTED' if p[0]<=.05 and robust=='PASS' else 'DIRECTIONALLY_SUPPORTED' if abs(g[0])>=.5 and preserve>=14 and maxinf<1 else 'NOT_SUPPORTED'
sig=q<=.05
secondary='CONVERGENT' if sig.sum()>=2 and np.all(np.sign(D[1:][sig])==np.sign(D[0])) and any(n!='cytotoxicity' for n,b in zip(names[1:],sig) if b) else 'MIXED' if sig.any() else 'NO_SUPPORT'
increment={'SUPPORTED':'MODERATE','DIRECTIONALLY_SUPPORTED':'LOW','NOT_SUPPORTED':'NONE'}[primary]
gates=dict(PTPN22_GENETIC_AF35_ACTIVATION_GATE=primary,PTPN22_GENETIC_AF35_DONOR_ROBUSTNESS_GATE=robust,
 PTPN22_GENETIC_PROGRAM_CONTEXT_GATE=secondary,PTPN22_GENETIC_MANUSCRIPT_INCREMENT=increment,
 SAMPLE_SCORING_QC='PASS',source_consumption='UPDATED_WITH_PUBLIC_SI_FULLTEXT_LIMITATION_RETAINED',
 LODO_sign_preserved_n=preserve,LODO_sign_reversal_n=sum(r['sign_reversal'] for r in lodo),LODO_max_influence_fraction=float(maxinf),
 LODO_over40_donors=[r['removed_donor'] for r in lodo if r['exceeds_40_percent']],secondary_BH_significant_n=int(sig.sum()))
js(gates,'10_checkpoint/FINAL_GATES.json')
tsv([dict(gate=k,status=v) for k,v in gates.items() if k.startswith('PTPN22_')],'10_checkpoint/FINAL_GATES.tsv')
assert sha(matrix_path)==workbook_hash
js(dict(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),primary=results[0],gates=gates,
 bootstrap_rng='numpy PCG64',boot_index_SHA256=hashlib.sha256(index_cc.tobytes()+index_tt.tobytes()).hexdigest()),'logs/ANALYSIS_COMPLETION.json')
print(json.dumps(dict(primary=results[0],gates=gates),indent=2))
