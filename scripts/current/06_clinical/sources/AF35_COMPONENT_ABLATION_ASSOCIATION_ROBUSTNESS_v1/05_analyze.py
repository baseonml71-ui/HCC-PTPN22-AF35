import importlib, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
f=importlib.import_module('01_freeze');OUT=f.OUT;P=f.P
for name in ['FREEZE_LOCK.json','SCORING_RECOVERY_LOCK.json']:
    for n,h in json.loads((OUT/name).read_text())['hashes'].items(): assert f.sha(OUT/n)==h
single=pd.read_csv(OUT/'AF35_SINGLE_GENE_SUBSET_MANIFEST.tsv',sep='\t')
random=pd.read_csv(OUT/'AF35_RANDOM_DROPOUT_SUBSET_MANIFEST.tsv',sep='\t')
genes=pd.read_csv(OUT/'AF35_ORDERED_GENES.tsv',sep='\t').gene.tolist()
def hedges(x,y):
    r=x[y==1];nr=x[y==0];n=len(x)
    sd=np.sqrt(((len(r)-1)*np.var(r,ddof=1)+(len(nr)-1)*np.var(nr,ddof=1))/(n-2))
    assert np.isfinite(sd) and sd>0
    return (1-3/(4*(n-2)-1))*(np.mean(r)-np.mean(nr))/sd
sg=[];rd=[];fullrows=[]
for c in ['ERP117672','GSE302495']:
    d=pd.read_csv(OUT/f'{c}_ANALYSIS_PATIENT_VALUES.tsv',sep='\t')
    y=d.response.to_numpy();full=d.AF35.to_numpy();gfull=hedges(full,y)
    old=pd.read_csv(P/f'05_results/tables/phase15c_public_benchmark/{c}_UNIVARIATE_BENCHMARK_EFFECTS.tsv',sep='\t')
    # The historical effect table uses metric as its identifier.
    key=next(k for k in ['metric','variable','candidate'] if k in old.columns)
    expected=old.loc[old[key].eq('AF35'),'hedges_g'].item()
    assert abs(gfull-expected)<1e-10, f'FULL_EFFECT_AUTHORITY_CONFLICT {c} {gfull} {expected}'
    fullrows.append(dict(cohort=c,n=len(d),responders=int(y.sum()),nonresponders=int((1-y).sum()),full_AF35_hedges_g=gfull,frozen_hedges_g=expected,abs_difference=abs(gfull-expected),status='PASS'))
    def calc(r):
        kept=r.retained_genes.split(';')
        assert kept==[g for g in genes if g in kept] and len(kept)==r.retained_gene_n
        score=d[kept].mean(axis=1).to_numpy()
        g=hedges(score,y)
        return dict(cohort=c,subset_id=r.subset_id,retained_gene_count=len(kept),retained_fraction=len(kept)/35,hedges_g=g,full_AF35_hedges_g=gfull,delta_g=g-gfull,absolute_delta_g=abs(g-gfull),direction='POSITIVE' if g>0 else ('NEGATIVE' if g<0 else 'ZERO'),score_rho=spearmanr(score,full).statistic,rho_population_n=len(d))
    for r in single.itertuples():
        row=calc(r);row['omitted_gene']=r.omitted_gene;row['gene_rank']=genes.index(r.omitted_gene)+1;sg.append(row)
    for r in random.loc[random.cohort.eq(c)].itertuples():
        row=calc(r);row.update(dropout_percent=r.dropout_percent,actual_removed_fraction=r.removed_gene_n/35,draw=r.draw);rd.append(row)
sg=pd.DataFrame(sg);rd=pd.DataFrame(rd)
sg.to_csv(OUT/'AF35_SINGLE_GENE_ASSOCIATION_ABLATION.tsv',sep='\t',index=False)
rd.to_csv(OUT/'AF35_RANDOM_DROPOUT_ASSOCIATION_ROBUSTNESS.tsv',sep='\t',index=False)
f.tsv(fullrows,'AF35_FULL_SCORE_EFFECT_REPRODUCTION.tsv')
ss=[];rs=[]
for c,s in sg.groupby('cohort',sort=False):
    ss.append(dict(cohort=c,positive_n=int((s.hedges_g>0).sum()),total_n=35,g_min=s.hedges_g.min(),g_max=s.hedges_g.max(),delta_g_min=s.delta_g.min(),delta_g_max=s.delta_g.max(),max_absolute_delta_g=s.absolute_delta_g.max(),min_score_rho=s.score_rho.min()))
for (c,pct),s in rd.groupby(['cohort','dropout_percent'],sort=False):
    q=np.quantile(s.hedges_g,[.025,.25,.5,.75,.975],method='linear')
    rs.append(dict(cohort=c,dropout_percent=pct,removed_gene_count=round(s.actual_removed_fraction.iloc[0]*35),retained_gene_count=int(s.retained_gene_count.iloc[0]),draw_n=len(s),g_p025=q[0],g_p25=q[1],median_hedges_g=q[2],g_p75=q[3],g_p975=q[4],median_delta_g=s.delta_g.median(),median_score_rho=s.score_rho.median(),positive_fraction=float((s.hedges_g>0).mean()),full_AF35_hedges_g=s.full_AF35_hedges_g.iloc[0],distribution_type='COMPONENT_PERTURBATION_NOT_PATIENT_CI'))
ss=pd.DataFrame(ss);rs=pd.DataFrame(rs)
ss.to_csv(OUT/'AF35_SINGLE_GENE_SUMMARY.tsv',sep='\t',index=False)
rs.to_csv(OUT/'AF35_RANDOM_DROPOUT_SUMMARY.tsv',sep='\t',index=False)
neg=sg.pivot(index='omitted_gene',columns='cohort',values='hedges_g')<0
f.tsv([dict(reverse_in_either_cohort=int(neg.any(axis=1).sum()),reverse_in_both_cohorts=int(neg.all(axis=1).sum()),zero_effect_deletions=int(sg.hedges_g.eq(0).sum()))],'AF35_CROSS_COHORT_DELETION_DIRECTIONS.tsv')
single_label='SINGLE_GENE_SENSITIVE' if (sg.hedges_g.le(0).any() or sg.absolute_delta_g.ge(.30).any()) else ('DISTRIBUTED' if sg.absolute_delta_g.le(.15).all() and sg.score_rho.ge(.95).all() else 'MODERATE_COMPONENT_SENSITIVITY')
primary=rs.loc[rs.dropout_percent.le(30)]
random_label='EARLY_COLLAPSE' if (primary.positive_fraction.lt(.80).any() or (primary.median_hedges_g<.5*primary.full_AF35_hedges_g).any()) else ('GRACEFUL_DEGRADATION' if rs.positive_fraction.ge(.95).all() and (rs.median_hedges_g>=.75*rs.full_AF35_hedges_g).all() and rs.median_score_rho.ge(.90).all() else 'INTERMEDIATE')
strong=ss.positive_n.eq(35).all() and primary.positive_fraction.ge(.95).all()
partial=any((ss.loc[ss.cohort.eq(c),'positive_n'].item()/35>=.80 and primary.loc[primary.cohort.eq(c),'positive_fraction'].ge(.80).all()) for c in ss.cohort)
cross='STRONG' if strong else ('PARTIAL' if partial else 'WEAK')
conditional=pd.read_csv(OUT/'AF35_PTPN22_CONDITIONAL_ASSOCIATION_SENSITIVITY.tsv',sep='\t')
cl='NOT_INTERPRETABLE' if not conditional.numerical_interpretability.eq('PASS').all() else ('POSITIVE_BOTH' if conditional.loc[conditional.term.eq('AF35'),'coefficient'].gt(0).all() else 'MIXED')
f.tsv([dict(domain=k,classification=v,basis='PROTOCOL_v1_PRE_OUTCOME_DESCRIPTIVE_RULE',claim_boundary='NO_OPTIMALITY_NO_CAUSAL_IMPORTANCE_NO_CLASSIFIER') for k,v in [('SINGLE_GENE_DEPENDENCE',single_label),('RANDOM_DROPOUT_BEHAVIOR',random_label),('CROSS_COHORT_DIRECTION_STABILITY',cross),('ANCHOR_READOUT_CONDITIONAL_SENSITIVITY',cl)]],'AF35_COMPONENT_ROBUSTNESS_INTERPRETATION.tsv')
print(ss.to_string(index=False));print(rs.to_string(index=False));print(conditional.to_string(index=False));print(single_label,random_label,cross,cl)
