from pathlib import Path
import importlib.util, json, hashlib, sys
import numpy as np
import pandas as pd
sys.dont_write_bytecode = True

OUT = Path(__file__).resolve().parents[1]
ROOT = OUT.parent
P = OUT / 'authority'
def write(df, name):
    pd.DataFrame(df).to_csv(OUT / name, sep='\t', index=False, na_rep='NA')
def read(rel):
    return pd.read_csv(P / rel, sep='\t')
def digest(p):
    with open(p,'rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()

spec = importlib.util.spec_from_file_location('frozen', P / '04_scripts/phase19/02_phase19a_patient_analysis.py')
frozen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(frozen)
rng = np.random.default_rng(20260908)
checks, strata, loo_rows, all_summaries = [], [], [], []

def check(name, ok, detail=''):
    checks.append(dict(check=name,status='PASS' if ok else 'FAIL',detail=detail))
    if not ok:
        write(checks,'ANALYSIS_VERIFICATION.tsv')
        raise AssertionError((name,detail))

def summary(v):
    return frozen.summarize_effects(dict(v), rng, draws=10000)

def estimate(df, flag, context, identity='AUTHOR_CLONE_ID'):
    patient_rows, summaries = [], []
    for score in ['PTPN22','AF35']:
        for patient,g in df.groupby('patient',sort=True):
            assert g[flag].notna().all()
            pos = g[flag].astype(bool)
            raw = g.loc[pos,score].mean()-g.loc[~pos,score].mean()
            pieces, supported = [], []
            for state,s in g.groupby('sub_cluster',sort=False):
                sp = s[flag].astype(bool)
                n1,n0 = int(sp.sum()),int((~sp).sum())
                d = s.loc[sp,score].mean()-s.loc[~sp,score].mean()
                w = min(n1,n0)
                strata.append(dict(identity=identity,context=context,score=score,patient=patient,state=state,
                                   positive_cells=n1,reference_cells=n0,weight=w,state_effect=d,
                                   included=bool(w),exclusion_reason='' if w else 'NO_WITHIN_STATE_COMPARATOR'))
                if w:
                    pieces.append((d,w)); supported.extend(s.index)
            adj = np.average([x[0] for x in pieces],weights=[x[1] for x in pieces]) if pieces else np.nan
            records = [{flag:bool(r[flag]),score:float(r[score]),'sub_cluster':r['sub_cluster']} for r in g.to_dict('records')]
            old = frozen.state_stratified_effect(records,flag,score)
            check(f'operator_{identity}_{context}_{score}_{patient}',np.isclose(adj,old,atol=1e-12,rtol=0,equal_nan=True))
            supported_df = g.loc[supported]
            sp = supported_df[flag].astype(bool)
            support_raw = supported_df.loc[sp,score].mean()-supported_df.loc[~sp,score].mean()
            patient_rows.append(dict(row_type='PATIENT',identity=identity,context=context,score=score,patient=patient,
                eligible_cells=len(g),positive_cells=int(pos.sum()),reference_cells=int((~pos).sum()),
                supported_states=len(pieces),eligible_states=g.sub_cluster.nunique(),supported_cells=len(supported),
                excluded_cells=len(g)-len(supported),raw_effect=raw,adjusted_effect=adj,
                supported_cell_raw_effect=support_raw,attenuation_difference=raw-adj,
                adjusted_to_raw_ratio=adj/raw if raw else np.nan,attenuation_fraction=1-adj/raw if raw else np.nan,
                biological_replicate='patient'))
        r = pd.DataFrame(patient_rows)
        r = r[r.score==score]
        valid = r.raw_effect.notna() & r.adjusted_effect.notna()
        for scale, col in [('RAW','raw_effect'),('STATE_ADJUSTED','adjusted_effect')]:
            effects = list(zip(r.loc[valid,'patient'],r.loc[valid,col]))
            s = dict(identity=identity,context=context,score=score,estimate=scale,**summary(effects))
            s['seed'] = 20260908
            s['biological_replicate']='patient'
            s['within_patient_weighting']='cells; state min-group-cell weights for adjusted'
            summaries.append(s); all_summaries.append(s)
            for patient,_ in effects:
                loo_rows.append(dict(identity=identity,context=context,score=score,estimate=scale,
                    omitted_patient=patient,remaining_patient_n=len(effects)-1,
                    mean_effect=np.mean([v for p,v in effects if p!=patient])))
        raw_mean,adj_mean = r.loc[valid,'raw_effect'].mean(),r.loc[valid,'adjusted_effect'].mean()
        for s in summaries[-2:]:
            s['raw_mean']=raw_mean; s['adjusted_mean']=adj_mean
            s['attenuation_difference']=raw_mean-adj_mean
            s['adjusted_to_raw_ratio']=adj_mean/raw_mean if raw_mean else np.nan
            s['attenuation_fraction']=1-adj_mean/raw_mean if raw_mean else np.nan
    return pd.DataFrame(patient_rows),pd.DataFrame(summaries)

def overlap_row(g,a,b,patient,tissue,pair,universe):
    x,y=g[a].astype(bool),g[b].astype(bool)
    n11=int((x&y).sum());n10=int((x&~y).sum());n01=int((~x&y).sum());n00=int((~x&~y).sum())
    union=n11+n10+n01
    return dict(patient=patient,tissue=tissue,pair=pair,context1=a,context2=b,eligibility=universe,
        status='JOINTLY_ELIGIBLE',eligible_clone_count=len(g),n11=n11,n10=n10,n01=n01,n00=n00,
        positive_union_n=union,jaccard=n11/union if union else np.nan,
        p_context2_given_context1=n11/(n11+n10) if n11+n10 else np.nan,
        p_context1_given_context2=n11/(n11+n01) if n11+n01 else np.nan,
        counting_unit='patient+tissue+author_clone_id',inference='DESCRIPTIVE_ONLY')

def main():
    replay= pd.read_csv(OUT / 'FROZEN_PHASE19_NUMERIC_REPRODUCTION.tsv',sep='\t')
    check('authority_replay',replay.status.eq('PASS').all())
    x=read('03_processed_data/phase3/GSE235863_BARCODE_CLONOTYPE_STATE_PTPN22.tsv.gz')
    # response_group is never used by new analyses.
    x=x.drop(columns='response_group')
    z=np.load(P/'03_processed_data/phase19a_af35_clonotype/GSE235863_AF35_RESPONSE_BLIND_PROJECTION.npz',allow_pickle=False)
    scores=pd.DataFrame({'barcode':z['barcode'],'AF35':z['af35_primary']})
    x=x.merge(scores,on='barcode',validate='one_to_one',how='left')
    x['PTPN22']=x.ptpn22_log1p_cp10k
    check('exact_link_58872',len(x)==58872 and x.barcode.nunique()==58872 and x.AF35.notna().all())
    check('ptpn22_normalization',np.allclose(x.PTPN22,np.log1p(10000*x.ptpn22_raw_count/x.total_raw_umi),atol=1e-12,rtol=0))
    check('20_released_states',x.sub_cluster.nunique()==20)
    check('author_size_exact',x.groupby('clonotype_key')['clone.size'].nunique().eq(1).all() and
          x.groupby('clonotype_key').size().eq(x.groupby('clonotype_key')['clone.size'].first()).all())
    x['expanded']=x['clone.size']>=2
    x['persistent']=x.persistent_within_tissue.map(lambda v: v is True or str(v).upper()=='TRUE')
    fate=read('05_results/tables/temporal/GSE235863_BASELINE_CLONE_FATE_TABLE.tsv')
    fate=fate[fate.analysis_role=='PRIMARY_BLOOD']
    check('fate_11411_unique',len(fate)==11411 and not fate.duplicated(['patient','clonotype_key']).any())
    fate=fate[['patient','clonotype_key','future_persistent_observed']].rename(columns={'future_persistent_observed':'later_observed'})
    x=x.merge(fate,on=['patient','clonotype_key'],how='left',validate='many_to_one')
    temporal=x[(x.timepoint=='pre')&(x.tissue=='P')&x.later_observed.notna()].copy()
    persistent=x[x.persistence_eligible].copy()
    check('n_patients_9_7_7',[x.patient.nunique(),persistent.patient.nunique(),temporal.patient.nunique()]==[9,7,7])
    universes={'EXPANSION':(x,'expanded'),'PERSISTENCE':(persistent,'persistent'),'LATER_OBSERVATION':(temporal,'later_observed')}
    patient_tables=[]
    for context,(g,flag) in universes.items():
        patients,s=estimate(g,flag,context)
        patient_tables.append(patients)
        write(pd.concat([patients,s.assign(row_type='SUMMARY',patient='AGGREGATE')],ignore_index=True),f'CLONE_STATE_ADJUSTED_{context}.tsv')
    primary=pd.concat(patient_tables,ignore_index=True)
    write(primary,'CLONE_STATE_PATIENT_EFFECTS.tsv')
    historic=read('05_results/tables/phase19a_af35_clonotype/AF35_CLONOTYPE_TEMPORAL_SENSITIVITY.tsv')
    rec=[]
    for context,analysis in [('EXPANSION','EXPANDED_VS_SINGLETON'),('PERSISTENCE','PERSISTENT_VS_NONPERSISTENT'),('LATER_OBSERVATION','BASELINE_FUTURE_OBSERVED')]:
        h=historic[(historic.analysis==analysis)&(historic.sensitivity_method=='AUTHOR_STATE_STRATIFIED')].iloc[0]
        vals=primary[(primary.context==context)&(primary.score=='AF35')].adjusted_effect
        check('af35_frozen_adjusted_'+context,np.isclose(vals.mean(),h.mean_patient_effect,atol=1e-12,rtol=0))
        rec.append(dict(context=context,mean_patient_effect=vals.mean(),frozen_mean=h.mean_patient_effect,
                        positive_n=int((vals>0).sum()),patient_n=len(vals),frozen_ci_low=h.bootstrap_ci_low,frozen_ci_high=h.bootstrap_ci_high))
    write(rec,'AF35_FROZEN_STATE_EFFECT_RECOVERY.tsv')
    ptp=read('05_results/tables/phase3/GSE235863_WITHIN_PATIENT_CLONAL_ASSOCIATIONS.tsv')
    for context,metric in [('EXPANSION','expanded_minus_singleton_ptpn22'),('PERSISTENCE','persistent_minus_nonpersistent_ptpn22')]:
        vals=primary[(primary.context==context)&(primary.score=='PTPN22')].raw_effect
        check('ptpn22_frozen_raw_'+context,np.isclose(vals.mean(),ptp[ptp.metric==metric].iloc[0].mean_patient_effect,atol=1e-12,rtol=0))
    # Clone-context rows preserve tissue-specific persistence eligibility.
    pc=persistent.groupby(['patient','tissue','clonotype_key'],as_index=False).agg(expanded=('expanded','first'),persistent=('persistent','first'))
    tc=temporal.groupby(['patient','tissue','clonotype_key'],as_index=False).agg(expanded=('expanded','first'),persistent=('persistent','first'),later_observed=('later_observed','first'))
    check('baseline_blood_persistence_joint_eligibility',temporal.persistence_eligible.all())
    check('temporal_persistence_identity',tc.persistent.eq(tc.later_observed).all())
    overlaps=[]
    for (pat,tissue),g in pc.groupby(['patient','tissue']):
        overlaps.append(overlap_row(g,'expanded','persistent',pat,tissue,'EXPANSION_PERSISTENCE','FROZEN_PERSISTENCE_ELIGIBLE_TISSUE_CLONES_ALL_TIMEPOINTS'))
    for pat,g in tc.groupby('patient'):
        overlaps.append(overlap_row(g,'expanded','later_observed',pat,'P','EXPANSION_LATER_OBSERVATION','BASELINE_PB_CLONES_WITH_EXPLICIT_FOLLOWUP'))
        overlaps.append(overlap_row(g,'persistent','later_observed',pat,'P','PERSISTENCE_LATER_OBSERVATION','BASELINE_PB_SHARED_ELIGIBILITY'))
    ov=pd.DataFrame(overlaps)
    check('overlap_counts',ov[['n11','n10','n01','n00']].sum(axis=1).eq(ov.eligible_clone_count).all())
    check('persistence_implies_expansion',not (pc.persistent&~pc.expanded).any())
    check('later_observed_implies_global_expansion',not (tc.later_observed.astype(bool)&~tc.expanded).any())
    write(ov,'CLONE_CONTEXT_ELIGIBILITY_AWARE_OVERLAP.tsv')
    write(pc,'CLONE_CONTEXT_PERSISTENCE_ELIGIBLE_CLONES.tsv')
    write(tc,'CLONE_CONTEXT_BASELINE_PB_SHARED_ELIGIBLE_CLONES.tsv')
    triples=[]
    for pat,g in tc.groupby('patient'):
        for e in [False,True]:
            for p in [False,True]:
                for y in [False,True]:
                    triples.append(dict(patient=pat,tissue='P',eligible_clone_n=len(g),expanded=e,persistent=p,later_observed=y,
                        clone_n=int(((g.expanded==e)&(g.persistent==p)&(g.later_observed==y)).sum()),
                        interpretation='SHARED_BASELINE_PB_ONLY; PERSISTENCE_EQUALS_LATER_OBSERVATION_BY_DEFINITION'))
    write(triples,'CLONE_CONTEXT_SHARED_UNIVERSE_TRIPLE_OVERLAP.tsv')
    # Patient + paired CDR3 identity: rederive size, persistence and later detection on identical cell universes.
    x['cdr3_key']=x.patient+'::'+x.cdr3_pair_aa
    check('cdr3_complete',x.cdr3_pair_aa.notna().all() and x.cdr3_pair_aa.str.split('|').map(len).eq(2).all())
    merge_audit=x.groupby(['patient','cdr3_pair_aa'],as_index=False).agg(author_clone_n=('clonotype_key','nunique'),cell_n=('barcode','size'))
    write(merge_audit[merge_audit.author_clone_n>1],'PAIRED_CDR3_AUTHOR_ID_MERGE_AUDIT.tsv')
    x['cdr3_expanded']=x.groupby('cdr3_key').barcode.transform('size')>=2
    pre=set(zip(x.loc[x.timepoint=='pre','patient'],x.loc[x.timepoint=='pre','tissue'],x.loc[x.timepoint=='pre','cdr3_pair_aa']))
    post=set(zip(x.loc[x.timepoint=='post','patient'],x.loc[x.timepoint=='post','tissue'],x.loc[x.timepoint=='post','cdr3_pair_aa']))
    x['cdr3_persistent']=[k in pre and k in post for k in zip(x.patient,x.tissue,x.cdr3_pair_aa)]
    x['cdr3_later']=[k in post for k in zip(x.patient,x.tissue,x.cdr3_pair_aa)]
    sens=[]; changed=[]
    for context,mask,flag,original in [
        ('EXPANSION',pd.Series(True,index=x.index),'cdr3_expanded','expanded'),
        ('PERSISTENCE',x.persistence_eligible,'cdr3_persistent','persistent'),
        ('LATER_OBSERVATION',(x.timepoint=='pre')&(x.tissue=='P')&x.later_observed.notna(),'cdr3_later','later_observed')]:
        g=x[mask].copy()
        for pat,r in g.groupby('patient'):
            changed.append(dict(context=context,patient=pat,eligible_cells=len(r),changed_flag_cells=int((r[flag]!=r[original]).sum())))
        patients,s=estimate(g,flag,context,'PAIRED_CDR3_ALPHA_BETA')
        sens.append(patients)
    write(pd.concat(sens,ignore_index=True),'PAIRED_CDR3_STATE_ADJUSTED_SENSITIVITY.tsv')
    write(changed,'PAIRED_CDR3_CONTEXT_CLASSIFICATION_AUDIT.tsv')
    # Matched compartment feasibility is computed before any expression contrast.
    meta=read('01_metadata/final_harmonization/GSE235863_TEMPORAL_COMPARTMENT_AUDIT.tsv')
    contexts=[]; clone_sizes=[]; co_cells=[]; source_rows=[]
    for (pat,time),g in x.groupby(['patient','timepoint']):
        if set(g.tissue)!={'P','T'}:
            contexts.append(dict(patient=pat,timepoint=time,status='UNMATCHED_COMPARTMENT',blood_cells=int((g.tissue=='P').sum()),tumor_cells=int((g.tissue=='T').sum())))
            continue
        for tissue,word in [('P','blood'),('T','liver tumor')]:
            sample=g.loc[g.tissue==tissue,'sample'].unique()
            check(f'single_sample_{pat}_{time}_{tissue}',len(sample)==1)
            m=meta[meta.sample_id==sample[0]]
            check(f'source_sample_{pat}_{time}_{tissue}',len(m)==1)
            m=m.iloc[0]
            title=f'{"Pre" if time=="pre" else "Post"}-treatment {word} of patient {pat} '
            check(f'source_title_{pat}_{time}_{tissue}',m.exact_geo_title.startswith(title))
            source_rows.append(m.to_dict())
        c=g.groupby('clonotype_key',as_index=False).agg(compartment_n=('tissue','nunique'),matched_pair_cell_n=('barcode','size'),global_cell_n=('clone.size','first'))
        c['co_detected']=c.compartment_n==2
        c['patient']=pat;c['timepoint']=time
        c['local_size_category']=pd.cut(c.matched_pair_cell_n,bins=[0,1,4,np.inf],labels=['1','2-4','>=5']).astype(str)
        clone_sizes.append(c)
        g=g.merge(c[['clonotype_key','co_detected','matched_pair_cell_n','local_size_category']],on='clonotype_key',validate='many_to_one')
        compartment_comparable=all(s.co_detected.nunique()==2 for _,s in g.groupby('tissue'))
        contexts.append(dict(patient=pat,timepoint=time,status='SOURCE_VALID_MATCHED',blood_cells=int((g.tissue=='P').sum()),tumor_cells=int((g.tissue=='T').sum()),
            eligible_clone_n=len(c),co_detected_clone_n=int(c.co_detected.sum()),restricted_clone_n=int((~c.co_detected).sum()),compartments_comparable=compartment_comparable,
            median_codetected_size=c.loc[c.co_detected,'matched_pair_cell_n'].median(),median_restricted_size=c.loc[~c.co_detected,'matched_pair_cell_n'].median(),
            timing_resolution='released_pre_post_phase; exact_same_day_not_verified'))
        if compartment_comparable: co_cells.append(g)
    ctx=pd.DataFrame(contexts)
    write(ctx,'TUMOR_BLOOD_MATCHED_CONTEXT_AUDIT.tsv')
    write(source_rows,'TUMOR_BLOOD_SOURCE_SAMPLE_BINDINGS.tsv')
    cs=pd.concat(clone_sizes,ignore_index=True)
    write(cs,'TUMOR_BLOOD_CODETECTION_ELIGIBLE_CLONES.tsv')
    size_audit=cs.groupby(['patient','timepoint','local_size_category','co_detected'],observed=True).agg(clone_n=('clonotype_key','size'),median_local_size=('matched_pair_cell_n','median'),median_global_size=('global_cell_n','median')).reset_index()
    write(size_audit,'TUMOR_BLOOD_CLONE_SIZE_CONFOUNDING_AUDIT.tsv')
    valid_patients=sorted({g.patient.iloc[0] for g in co_cells})
    feasible=len(valid_patients)>=3
    (OUT/'TUMOR_BLOOD_CODETECTION_FEASIBILITY.md').write_text(f'''# 肿瘤/外周血共检出可行性

TUMOR_BLOOD_CODETECTION_FEASIBILITY_GATE = {'PASS' if feasible else 'INSUFFICIENT'}

来源有效同患者同 pre/post 阶段匹配 {int((ctx.status=='SOURCE_VALID_MATCHED').sum())} 个上下文，表达比较可用 {len(valid_patients)} 位患者。每个组织一个明确释放样本，GEO 标题逐行核验，详见 SOURCE_SAMPLE_BINDINGS。只验证同治疗阶段，不声称同一天采样。保留 P27 GEO characteristics 患者字段冲突；沿用冻结标题/样本 token 映射。
匹配上下文内共有 {len(cs)} 个克隆上下文记录，共检出 {int(cs.co_detected.sum())}，局限检出 {int((~cs.co_detected).sum())}；跨时点不可当成独立克隆或独立患者。
共检出至少需要两个观测细胞，所以 size=1 必然局限检出；克隆大小分布和分层计数见 CLONE_SIZE_CONFOUNDING_AUDIT。测序深度和组织可采样性仍可影响检测，不据此作交通、迁移或播散推断。
技术 PASS 的标准已在结果计算前写入方案，不要求任何表达方向或显著性。通过后分别计算患者×时点×组织内两组细胞均值差及同一状态校正差，组织等权后时点等权，最后患者等权 bootstrap；Supplement-only。
''',encoding='utf-8')
    co_patient=[]; co_detail=[]
    if feasible:
        for score in ['PTPN22','AF35']:
            for g in co_cells:
                for tissue,t in g.groupby('tissue'):
                    pos=t.co_detected.astype(bool)
                    records=[{'co_detected':bool(r.co_detected),score:float(getattr(r,score)),'sub_cluster':r.sub_cluster} for r in t.itertuples()]
                    co_detail.append(dict(patient=g.patient.iloc[0],timepoint=g.timepoint.iloc[0],tissue=tissue,score=score,
                        positive_cells=int(pos.sum()),reference_cells=int((~pos).sum()),
                        raw_effect=t.loc[pos,score].mean()-t.loc[~pos,score].mean(),
                        adjusted_effect=frozen.state_stratified_effect(records,'co_detected',score)))
        detail=pd.DataFrame(co_detail)
        check('codetection_complete_state_comparators',detail.adjusted_effect.notna().all())
        pair=detail.groupby(['patient','timepoint','score'],as_index=False)[['raw_effect','adjusted_effect']].mean()
        patient=pair.groupby(['patient','score'],as_index=False)[['raw_effect','adjusted_effect']].mean()
        patient['matched_timepoints']=patient.patient.map(pair.groupby('patient').timepoint.nunique())
        patient['weighting']='equal_compartment_then_equal_timepoint_then_equal_patient'
        write(detail,'TUMOR_BLOOD_CODETECTION_COMPARTMENT_EFFECTS.tsv')
        write(pair,'TUMOR_BLOOD_CODETECTION_MATCHED_PAIR_EFFECTS.tsv')
        write(patient,'TUMOR_BLOOD_CODETECTION_PATIENT_EFFECTS.tsv')
        for score,g in patient.groupby('score'):
            for label,col in [('RAW','raw_effect'),('STATE_ADJUSTED','adjusted_effect')]:
                all_summaries.append(dict(identity='AUTHOR_CLONE_ID',context='TUMOR_BLOOD_CODETECTION',score=score,estimate=label,
                    **summary(zip(g.patient,g[col])),seed=20260908,biological_replicate='patient',within_patient_weighting='equal_compartment_then_equal_timepoint'))
    write(all_summaries,'CLONE_STATE_EFFECT_SUMMARY.tsv')
    write(strata,'CLONE_STATE_STRATUM_EFFECTS_AND_SUPPORT.tsv')
    write(loo_rows,'CLONE_STATE_PATIENT_OMISSION.tsv')
    write(checks,'ANALYSIS_VERIFICATION.tsv')
    versions={'python':sys.version,'executable':sys.executable,'numpy':np.__version__,'pandas':pd.__version__,'seed':20260908,'bootstrap_draws':10000}
    (OUT/'PYTHON_SESSION_INFO.json').write_text(json.dumps(versions,indent=2),encoding='utf-8')
    (OUT/'CLONE_STATE_ADJUSTMENT_AUTHORITY_AUDIT.md').write_text('''# 状态校正权威实现审计

CLONE_STATE_ADJUSTMENT_AUTHORITY_GATE = PASS

权威实现：HCC_ICI_project/04_scripts/phase19/02_phase19a_patient_analysis.py，state_stratified_effect()、group_effect()、summarize_effects() 及 main() 内的集合定义。
58,872 个 barcode 一一链接；AF35 使用冻结 35 基因 gene-wise z 等权均值，不含 PTPN22。本轮不重新计算或训练 AF35。
基线细胞集合：timepoint=pre、tissue=P，且在 analysis_role=PRIMARY_BLOOD 的冻结 fate 表内；7 位患者。状态字段是 TCR 作者 sub_cluster 原始值（20 个），不是 H5AD 别名、重新聚类或合并状态。
每个患者每个状态计算 positive minus reference 的细胞均值差 d_ps，仅保留两组均存在的层；患者校正差=sum(min(n1,n0)*d_ps)/sum(min(n1,n0))。没有回归、残差化、患者内重新 z 标准化或细胞 bootstrap。患者间等权均值和中位数；10,000 次患者有放回 bootstrap 均值的 2.5/97.5 百分位，逐患者遗漏均值。
原始完整脚本在本目录 reproduction 内重放，输出路径全部重定向；种子 20260831、原 RNG 调用顺序及全部敏感性调用保留。586 项数值结果与历史 JSON 在绝对误差 1e-12 内相同。
AF35 时间状态校正精确结果：6/7 正向，mean=0.031249209370383873，CI=[-0.0052920651712841955,0.057140672748183234]。扩增和持续检出 AF35 状态分层汇总已在 v19 存在，本轮只恢复、细化并核对，不宣称新发现。
本轮新增区间用相同患者 bootstrap 算法、预先固定 seed=20260908，不沿用旧随机流的任意中间位置；旧 CI 另表原样保留。全部本轮患者状态差另以 pandas 分层计算，逐项与原函数交叉核验。
原脚本没有状态内两组共存的层会被排除。当前输出公开层权重、被排除细胞数及共同状态支持集的未分层差，因此 raw-adjusted 是权重/状态组成相关的描述差，不是中介比例或因果分解。
PTPN22 v6 时间主结果的克隆等权及患者中位数不等于本轮细胞加权+患者均值；本轮时间面板标记新同算子敏感性，不替换历史结果。
历史文件 SHA256 绑定见 HISTORICAL_SOURCE_HASH_VERIFICATION.tsv；投影/映射哈希与冻结 contract 相符。protected-input 审计覆盖指定输入和正文/Figure 3 文件，不声称锁定整个工作区。
''',encoding='utf-8')
    (OUT/'ANALYSIS_STAGE_STATUS.json').write_text(json.dumps({'authority':'PASS','expansion':'PASS','persistence':'PASS','overlap':'PASS',
        'codetection_feasibility':'PASS' if feasible else 'INSUFFICIENT','codetection_analysis':'PASS' if feasible else 'NOT_RUN',
        'analysis_checks':len(checks),'codetection_patients':len(valid_patients)},indent=2),encoding='utf-8')
    print(pd.DataFrame(all_summaries)[['identity','context','score','estimate','eligible_patient_n','mean_patient_effect','bootstrap_ci_low','bootstrap_ci_high','positive_patient_n']].to_string(index=False))

if __name__=='__main__': main()
