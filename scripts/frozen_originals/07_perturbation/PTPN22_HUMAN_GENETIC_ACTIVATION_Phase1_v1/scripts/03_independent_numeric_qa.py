from pathlib import Path
import json,math,statistics,itertools,hashlib
import numpy as np
import pandas as pd
from scipy.stats import permutation_test,false_discovery_control
ROOT=Path(__file__).resolve().parents[1];I=ROOT/'01_manifest/inputs'
checks=[]
def check(name,value,limit=1e-10):
    ok=bool(value<=limit);checks.append(dict(check=name,max_abs_error=float(value),tolerance=limit,status='PASS' if ok else 'FAIL'));assert ok,name
# 独立标量评分：statistics.mean/stdev，不调用主分析score_fixed。
m=pd.read_excel(I/'GSE94859_RPKMs.xlsx',engine='openpyxl');cols=m.columns[2:].tolist()
defs=pd.read_csv(I/'PTPN22_MECHANISM_PROGRAMS_FROZEN.tsv',sep='\t')
export=pd.read_csv(ROOT/'02_scores/FROZEN_SCORES_BEFORE_GENOTYPE_JOIN.tsv',sep='\t').set_index('processed_matrix_column').loc[cols]
for program,sub in defs.groupby('program',sort=False):
    name='AF35' if program=='PTPN22_activation_feedback_35' else program
    rows=[]
    for gene in sub.gene:
        vals=m.loc[m.Gene_symbol==gene].iloc[0,2:].to_numpy(float)
        logvals=[math.log2(v+1) for v in vals];mu=statistics.mean(logvals);sd=statistics.stdev(logvals)
        rows.append([(v-mu)/sd for v in logvals])
    mean_score=[statistics.mean(row[j] for row in rows) for j in range(32)]
    mu=statistics.mean(mean_score);sd=statistics.stdev(mean_score)
    independent=np.array([(v-mu)/sd for v in mean_score])
    check(f'independent_scalar_scoring_{name}',np.max(np.abs(independent-export[name].to_numpy())))
primary=pd.read_csv(ROOT/'03_primary/GSE94859_AF35_GENOTYPE_INTERACTION.tsv',sep='\t').iloc[0]
d=pd.read_csv(ROOT/'03_primary/GSE94859_AF35_ACTIVATION_DELTA.tsv',sep='\t')
cc=d.loc[d.genotype=='CC','Delta_AF35'].to_numpy();tt=d.loc[d.genotype=='TT','Delta_AF35'].to_numpy()
check('primary_difference_from_donor_export',abs(statistics.mean(tt)-statistics.mean(cc)-primary.D_TT_minus_CC))
res=permutation_test((tt,cc),lambda a,b:np.mean(a)-np.mean(b),n_resamples=np.inf,alternative='two-sided',permutation_type='independent',vectorized=False)
assert len(res.null_distribution)==12870
check('independent_scipy_exact_permutation',abs(res.pvalue-primary.exact_P))
J=math.gamma(7)/(math.sqrt(7)*math.gamma(6.5));sp=math.sqrt((7*statistics.variance(cc)+7*statistics.variance(tt))/14)
check('independent_Hedges_g',abs(J*(statistics.mean(tt)-statistics.mean(cc))/sp-primary.Hedges_g))
sec=pd.read_csv(ROOT/'04_secondary/GSE94859_PROGRAM_GENOTYPE_INTERACTION.tsv',sep='\t')
check('independent_scipy_BH',np.max(np.abs(false_discovery_control(sec.exact_P.to_numpy(),method='bh')-sec.BH_q)))
boot=pd.read_csv(ROOT/'08_qa/AF35_BOOTSTRAP_DISTRIBUTION.tsv',sep='\t')
check('bootstrap_D_CI_from_saved_distribution',max(abs(boot.D_AF35.quantile(.025)-primary.bootstrap_D_CI95_low),abs(boot.D_AF35.quantile(.975)-primary.bootstrap_D_CI95_high)))
check('bootstrap_g_CI_from_saved_distribution',max(abs(boot.Hedges_g_AF35.quantile(.025)-primary.bootstrap_g_CI95_low),abs(boot.Hedges_g_AF35.quantile(.975)-primary.bootstrap_g_CI95_high)))
loo=pd.read_csv(ROOT/'05_robustness/GSE94859_AF35_LODO.tsv',sep='\t');err=[]
for r in loo.itertuples():
    s=d[d.individual_id!=r.removed_donor]
    estimate=s.loc[s.genotype=='TT','Delta_AF35'].mean()-s.loc[s.genotype=='CC','Delta_AF35'].mean()
    err.append(abs(estimate-r.LODO_D_AF35))
check('independent_LODO_from_donor_export',max(err))
meta=pd.read_csv(I/'GSE94859_SAMPLE_PAIRING_LEDGER.tsv',sep='\t').drop_duplicates('individual_id')
s=d.merge(meta[['individual_id','age','sex']],on='individual_id',validate='one_to_one')
X=np.column_stack([np.ones(16),(s.genotype=='TT').astype(int),s.age-s.age.mean(),(s.sex=='Male').astype(int)])
beta=np.linalg.lstsq(X,s.Delta_AF35,rcond=None)[0]
sen=pd.read_csv(ROOT/'05_robustness/GSE94859_AF35_COVARIATE_SENSITIVITY.tsv',sep='\t')
check('sensitivity_QR_lstsq_coefficient',np.max(np.abs(beta-sen.beta.to_numpy())))
checks.append(dict(check='no_p_zero_all_nine_endpoints',max_abs_error=0,tolerance=0,status='PASS'))
assert primary.exact_P>0 and (sec.exact_P>0).all()
pd.DataFrame(checks).to_csv(ROOT/'08_qa/INDEPENDENT_NUMERIC_QA.tsv',sep='\t',index=False)
print('Independent numeric QA passed:',len(checks),'checks; scipy exact P:',res.pvalue)
