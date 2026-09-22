from pathlib import Path
import json,hashlib,datetime
import numpy as np,pandas as pd
from coupling_math import decompose,technical_residual
ROOT=Path(__file__).resolve().parents[1];OLD=ROOT.parent/'PTPN22_AF35_COUPLING_ARCHITECTURE_v1'
def read(p):return pd.read_csv(p,sep='\t')
def save(p,rows):pd.DataFrame(rows).to_csv(ROOT/p,sep='\t',index=False,na_rep='NA',float_format='%.15g')
lock=json.loads((ROOT/'02_manifest/FREEZE_LOCK.json').read_text())
for p,h in lock['hashes'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
assert not (ROOT/'11_checkpoint/FINAL_GATES.json').exists()
(ROOT/'logs/EFFECT_RUN_STARTED.json').write_text(json.dumps(dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),manifest_sha256=lock['hashes']['02_manifest/PTPN22_STATE_AXIS_SPECIFICITY_ANALYSIS_MANIFEST_V1.md']),indent=2))
data=read(ROOT/'data/ANALYSIS_INPUT_RESPONSE_BLIND.tsv.gz');mapping=read(ROOT/'03_lineage/PTPN22_STATE_LINEAGE_MAP.tsv');states=mapping.released_state.tolist();smap=mapping.set_index('released_state').lineage.to_dict();data['lineage']=data.sub_cluster.map(smap)
patients=sorted(data.patient.unique(),key=lambda s:int(s[1:]));identities=[];lopo=[]
S='STATE_SEPARATION_PREDOMINANT';W='WITHIN_STATE_PREDOMINANT';U='COUPLING_ARCHITECTURE_UNRESOLVED';M='MIXED_COUPLING_ARCHITECTURE';NI='NOT_INTERPRETABLE'
def sign(x):return int(np.sign(x)) if np.isfinite(x) and abs(x)>1e-12 else 0
def one(f,patient,analysis,lineage='ALL',labels=None,adjust=False):
    x=f.ptpn22_log1p_cp10k.to_numpy();y=f.AF35.to_numpy();labs=f.sub_cluster.to_numpy() if labels is None else labels
    row=dict(source='GSE235863',patient=patient,analysis=analysis,lineage=lineage,cell_n=len(f),state_n=len(np.unique(labs)),PTPN22_variance=float(np.var(x)) if len(x) else np.nan,AF35_variance=float(np.var(y)) if len(y) else np.nan)
    reasons=[]
    if len(f)<50:reasons.append('LT50_CELLS')
    if row['state_n']<2:reasons.append('LT2_STATES')
    if not np.isfinite(x).all() or not np.isfinite(y).all():reasons.append('NONFINITE')
    if row['PTPN22_variance']==0:reasons.append('PTPN22_ZERO_VARIANCE')
    if row['AF35_variance']==0:reasons.append('AF35_ZERO_VARIANCE')
    if adjust and not reasons:x,y,rank=technical_residual(x,y,f.total_UMI,f.detected_genes);row['design_rank']=rank
    out=None if reasons else decompose(x,y,labs)
    if out is None:
        if not reasons:reasons.append('ADJUSTED_ZERO_VARIANCE')
        row.update(eligible=False,reason=';'.join(reasons))
        row.update({k:np.nan for k in ['C_total','C_between','C_within','identity_error','residual_Pearson','residual_Spearman','between_abs_gt_within','sign_total','sign_between','sign_within']});row['geometry']='NOT_EVALUABLE';return row,[]
    detail=out.pop('states');row.update(out);row.update(eligible=True,reason='ELIGIBLE');t,b,w=[out[k] for k in ['C_total','C_between','C_within']]
    row.update(between_abs_gt_within=abs(b)>abs(w),sign_total=sign(t),sign_between=sign(b),sign_within=sign(w),geometry='TOTAL_NEAR_ZERO' if abs(t)<.02 else 'OPPOSING_COMPONENTS' if b*w<0 else 'CONCORDANT_POSITIVE' if t>0 else 'CONCORDANT_NEGATIVE')
    identities.append({k:row[k] for k in ['source','patient','analysis','lineage','cell_n','C_total','C_between','C_within','identity_error']}|dict(PASS=abs(out['identity_error'])<1e-12))
    return row,detail
def core(f,min_n=3):
    f=f[f.eligible];n=len(f)
    if n<min_n:return U
    sg=sign(f.C_total.median());b=f.C_between.abs().median();w=f.C_within.abs().median()
    if abs(f.C_total.median())<.02 or (f.C_total.map(sign)==sg).mean()<2/3:return U
    if b>=.02 and w>=.02 and (f.C_between*f.C_within<0).mean()>=2/3:return 'OPPOSING_COUPLING_ARCHITECTURE'
    if (f.C_between.abs()>f.C_within.abs()).mean()>=2/3 and b>=1.5*w and (f.residual_Pearson.abs().median()<.10 or w<=b/1.5):return S
    if (f.C_within.abs()>f.C_between.abs()).mean()>=2/3 and w>=1.5*b and (f.residual_Pearson.map(sign)==sg).mean()>=2/3 and (f.residual_Spearman.map(sign)==sg).mean()>=2/3:return W
    return M if b>=.02 and w>=.02 else U
def summary(rows,label,lineage='ALL'):
    f=pd.DataFrame(rows);f=f[f.eligible];n=len(f);cg=core(f,4);stable=n>=4
    for p in f.patient:
        z=f[f.patient!=p];lg=core(z);ok=lg==cg and sign(z.C_total.median())==sign(f.C_total.median());stable=stable and ok
        lopo.append(dict(analysis=label,lineage=lineage,left_out_patient=p,eligible_patient_n=len(z),full_core_gate=cg,lopo_core_gate=lg,qualitative_stable=ok,median_C_total=z.C_total.median(),median_C_between=z.C_between.median(),median_C_within=z.C_within.median()))
    final=NI if n<4 else cg if stable else U
    lineagegate=NI if final==NI else 'MIXED' if final==M else 'UNRESOLVED' if final in [U,'OPPOSING_COUPLING_ARCHITECTURE'] else final
    b=f.C_between.abs().median();w=f.C_within.abs().median()
    return dict(source='GSE235863',analysis=label,lineage=lineage,eligible_patient_n=n,positive_total_n=int((f.C_total>1e-12).sum()),between_dominant_n=int((f.C_between.abs()>f.C_within.abs()).sum()),opposing_n=int((f.C_between*f.C_within<0).sum()),median_C_total=f.C_total.median(),median_C_between=f.C_between.median(),median_C_within=f.C_within.median(),median_abs_between=b,median_abs_within=w,median_abs_total=f.C_total.abs().median(),median_abs_between_over_within=b/w if w>0 else np.nan,median_residual_Pearson=f.residual_Pearson.median(),median_residual_Spearman=f.residual_Spearman.median(),median_abs_residual_Pearson=f.residual_Pearson.abs().median(),core_gate=cg,final_gate=final,lineage_gate=lineagegate,LOPO_stable=stable)

full=[];contributions=[];lineageonly=[]
prior=read(OLD/'03_primary/PTPN22_AF35_PATIENT_COUPLING_DECOMPOSITION.tsv').set_index('patient')
for patient in patients:
    f=data[data.patient==patient];r,detail=one(f,patient,'FULL_20STATE');full.append(r)
    for k in ['C_total','C_between','C_within']:assert abs(r[k]-prior.loc[patient,k])<1e-10
    lookup={q['state']:q for q in detail}
    for state in states:
        d=lookup.get(state);v=d['between_contribution'] if d else 0.
        contributions.append(dict(source='GSE235863',patient=patient,state=state,lineage=smap[state],cell_n=d['cell_n'] if d else 0,cell_fraction=d['weight'] if d else 0,mean_Xz=d['mean_Xz'] if d else np.nan,mean_Yz=d['mean_Yz'] if d else np.nan,signed_contribution=v,absolute_contribution=abs(v),full_patient_eligible=r['eligible']))
    lr,_=one(f,patient,'LINEAGE_ONLY',labels=f.lineage.to_numpy());lr['C_between_20state']=r['C_between'];lr['C_between_lineage']=lr['C_between'];lr['C_within_lineage']=lr['C_within'];lr['finer_state_within_lineage']=r['C_between']-lr['C_between']
    safe=abs(lr['C_total'])>=.02 and sign(lr['C_between'])==sign(lr['C_within'])==sign(lr['C_total'])!=0
    lr['lineage_total_share']=lr['C_between']/lr['C_total'] if safe else np.nan
    safe=abs(lr['C_total'])>=.02 and abs(r['C_between'])>=.02 and sign(lr['C_between'])==sign(lr['finer_state_within_lineage'])==sign(r['C_between'])!=0
    lr['lineage_share_of_B20']=lr['C_between']/r['C_between'] if safe else np.nan
    lineageonly.append(lr)
baseline=summary(full,'FULL_20STATE');assert baseline['final_gate']==S
primary=[];technical=[];rare=[];eligibility=[];summaries=[]
for lineage in ['CD4','CD8']:
    a=[];b=[];c=[]
    for patient in patients:
        f=data[(data.patient==patient)&(data.lineage==lineage)]
        r,_=one(f,patient,'PRIMARY',lineage);a.append(r);primary.append(r)
        eligibility.append({k:r[k] for k in ['source','patient','lineage','cell_n','state_n','PTPN22_variance','AF35_variance','eligible','reason']})
        tr,_=one(f,patient,'TECHNICAL_ADJUSTED',lineage,adjust=True);technical.append(tr);b.append(tr)
        kept=f[f.sub_cluster.map(f.sub_cluster.value_counts())>=20];rr,_=one(kept,patient,'RARE_STATE_GE20',lineage);rr['original_cell_n']=len(f);rr['retained_fraction']=len(kept)/len(f);rare.append(rr);c.append(rr)
    summaries.extend([summary(a,'PRIMARY',lineage),summary(b,'TECHNICAL_ADJUSTED',lineage),summary(c,'RARE_STATE_GE20',lineage)])
    save(f'03_lineage/PTPN22_AF35_{lineage}_COUPLING_DECOMPOSITION.tsv',pd.DataFrame(a).assign(**{f'C_total_{lineage}':[r['C_total'] for r in a],f'C_between_{lineage}':[r['C_between'] for r in a],f'C_within_{lineage}':[r['C_within'] for r in a]}))
loso=[];lososummary=[]
for removed in states+['REMOVE_ALL_CD4','REMOVE_ALL_CD8']:
    removed_states=[s for s in states if smap[s]==removed.replace('REMOVE_ALL_','')] if removed.startswith('REMOVE_ALL_') else [removed]
    rows=[]
    for patient in patients:
        orig=data[data.patient==patient];f=orig[~orig.sub_cluster.isin(removed_states)]
        row,_=one(f,patient,'LOSO_'+removed);row.update(removed_state=removed,removed_cell_n=len(orig)-len(f));rows.append(row);loso.append(row)
    s=summary(rows,'LOSO_'+removed);s['removed_state']=removed;s['removal_type']='LINEAGE' if len(removed_states)>1 else 'SINGLE_STATE'
    for metric in ['C_total','C_between','C_within','abs_between']:
        key='median_'+metric;old=baseline[key];s['full_'+key]=old;s['signed_change_'+metric]=s[key]-old;s['absolute_change_'+metric]=abs(s[key]-old);s['relative_change_'+metric]=(s[key]-old)/abs(old) if abs(old)>=.02 else np.nan
    eligible=[r['patient'] for r in rows if r['eligible']];bf=pd.DataFrame(full).set_index('patient').loc[eligible]
    s['matched_reference_patient_n']=len(bf);s['matched_reference_median_abs_between']=bf.C_between.abs().median()
    s['between_reduction_fraction']=(baseline['median_abs_between']-s['median_abs_between'])/baseline['median_abs_between'];s['architecture_changed']=s['final_gate']!=S
    s['HIGH_INFLUENCE_STATE']=s['architecture_changed'] or s['between_reduction_fraction']>=.30
    s['catastrophic_collapse']=s['between_reduction_fraction']>=.50 or s['median_abs_between']<.02 or sign(s['median_C_between'])!=sign(baseline['median_C_between']) or s['eligible_patient_n']<4
    lososummary.append(s)
save('03_lineage/PTPN22_LINEAGE_ELIGIBILITY.tsv',eligibility)
save('03_lineage/PTPN22_AF35_LINEAGE_RESIDUAL_COUPLING.tsv',[{k:r[k] for k in ['source','patient','lineage','cell_n','state_n','eligible','reason','residual_Pearson','residual_Spearman']} for r in primary])
save('03_lineage/LINEAGE_GATE_SUMMARY.tsv',summaries)
save('04_lineage_only/PTPN22_AF35_LINEAGE_ONLY_DECOMPOSITION.tsv',lineageonly)
save('05_loso/PTPN22_AF35_LOSO_BY_PATIENT.tsv',loso)
save('05_loso/PTPN22_AF35_LOSO_SUMMARY.tsv',lososummary)
save('06_contributions/PTPN22_AF35_STATE_BETWEEN_CONTRIBUTION.tsv',contributions)
cs=[];cdf=pd.DataFrame(contributions)
for state in states:
    z=cdf[(cdf.state==state)&cdf.full_patient_eligible];cs.append(dict(source='GSE235863',state=state,lineage=smap[state],eligible_patient_n=len(z),observed_patient_n=int((z.cell_n>0).sum()),median_signed_contribution=z.signed_contribution.median(),median_absolute_contribution=z.absolute_contribution.median(),min_signed_contribution=z.signed_contribution.min(),max_signed_contribution=z.signed_contribution.max(),min_absolute_contribution=z.absolute_contribution.min(),max_absolute_contribution=z.absolute_contribution.max(),positive_n=int((z.signed_contribution>1e-12).sum()),negative_n=int((z.signed_contribution< -1e-12).sum()),zero_n=int((z.signed_contribution.abs()<=1e-12).sum())))
save('06_contributions/PTPN22_AF35_STATE_BETWEEN_CONTRIBUTION_SUMMARY.tsv',cs)
save('07_sensitivity/PTPN22_AF35_LINEAGE_TECHNICAL_ADJUSTMENT.tsv',technical)
save('07_sensitivity/PTPN22_AF35_LINEAGE_RARE_STATE_SENSITIVITY.tsv',rare)
save('07_sensitivity/PATIENT_LOPO_ALL_ANALYSES.tsv',lopo)
save('09_qa/DECOMPOSITION_IDENTITY_AUDIT.tsv',identities)
save('09_qa/FULL_20STATE_REPRODUCTION.tsv',full)
lookup={(r['lineage'],r['analysis']):r for r in summaries}
preserved=[];supported=[];attenuation={}
for lin in ['CD4','CD8']:
    a=lookup[lin,'PRIMARY'];b=lookup[lin,'TECHNICAL_ADJUSTED'];c=lookup[lin,'RARE_STATE_GE20']
    ok=a['lineage_gate'] in [S,W,'MIXED'] and a['lineage_gate']==b['lineage_gate'] and sign(a['median_C_total'])==sign(b['median_C_total'])
    if ok:preserved.append(lin)
    if a['lineage_gate']==b['lineage_gate']==c['lineage_gate']==S and a['median_C_total']>0 and b['median_C_total']>0:supported.append(lin)
    attenuation[lin]=1-b['median_abs_total']/a['median_abs_total'] if a['median_abs_total']>0 else None
tg='ARCHITECTURE_PRESERVED' if len(preserved)==2 else 'PARTIALLY_PRESERVED' if preserved else 'NOT_PRESERVED'
single=pd.DataFrame(lososummary);single=single[single.removal_type=='SINGLE_STATE']
lg='FAIL' if single.catastrophic_collapse.any() or single.architecture_changed.sum()>=10 else 'PARTIAL' if single.architecture_changed.any() else 'PASS'
ld=pd.DataFrame(lineageonly);bratio=ld.C_between_lineage.abs().median()/baseline['median_abs_between'];fine=ld.finer_state_within_lineage.abs().median();agree=(ld.C_between_lineage.map(sign)==ld.C_between_20state.map(sign)).mean()
material=ld.C_between_lineage.abs().median()>=.02 and agree>=2/3 and bratio>=.25
dominant=material and bratio>=.75 and fine<=.25*baseline['median_abs_between']
ag=[lookup[lin,'PRIMARY'] for lin in ['CD4','CD8']]
if any(a['eligible_patient_n']<4 for a in ag):branch='STATE_AXIS_SPECIFICITY_UNRESOLVED'
elif all(a['lineage_gate']!=S for a in ag) and dominant:branch='LINEAGE_SEPARATION_DOMINANT'
elif supported and lg=='PASS':branch='MIXED_LINEAGE_AND_STATE_ARCHITECTURE' if material else 'WITHIN_LINEAGE_STATE_ARCHITECTURE_SUPPORTED'
else:branch='STATE_AXIS_SPECIFICITY_UNRESOLVED'
increment='NONE' if branch in ['LINEAGE_SEPARATION_DOMINANT','STATE_AXIS_SPECIFICITY_UNRESOLVED'] else 'MODERATE'
strong=False
for lin in supported:
    a=lookup[lin,'PRIMARY'];other=lookup['CD8' if lin=='CD4' else 'CD4','PRIMARY']
    strong=strong or (a['eligible_patient_n']>=5 and a['between_dominant_n']/a['eligible_patient_n']>=.8 and a['median_abs_between_over_within']>=3 and a['median_abs_between']>=.10 and other['median_C_total']>=0 and other['lineage_gate']!=W)
if increment=='MODERATE' and lg=='PASS' and tg=='ARCHITECTURE_PRESERVED' and all(lookup[lin,'PRIMARY']['eligible_patient_n']>=5 and attenuation[lin]<.30 for lin in supported) and (len(supported)==2 or strong):increment='HIGH'
gates=dict(PTPN22_CD4_WITHIN_LINEAGE_GATE=lookup['CD4','PRIMARY']['lineage_gate'],PTPN22_CD8_WITHIN_LINEAGE_GATE=lookup['CD8','PRIMARY']['lineage_gate'],PTPN22_LOSO_ROBUSTNESS_GATE=lg,PTPN22_TECHNICAL_ADJUSTMENT_GATE=tg,PTPN22_STATE_AXIS_SPECIFICITY_GATE=branch,STATE_AXIS_MANUSCRIPT_INCREMENT=increment,supported_lineages=supported,technical_preserved_lineages=preserved,technical_total_attenuation=attenuation,lineage_between_magnitude_ratio=bratio,lineage_between_median=float(ld.C_between_lineage.median()),fine_state_within_lineage_median_abs=float(fine),lineage_material=bool(material),lineage_dominant=bool(dominant),single_state_architecture_changes=int(single.architecture_changed.sum()),single_state_high_influence=single.loc[single.HIGH_INFLUENCE_STATE,'removed_state'].tolist(),single_state_catastrophic=single.loc[single.catastrophic_collapse,'removed_state'].tolist(),max_between_reduction_fraction=float(single.between_reduction_fraction.max()),patient_n=6,cell_n=9088,independent_replication='NOT_INTERPRETABLE',outcome_used=False,causal_claim=False,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(ROOT/'11_checkpoint/FINAL_GATES.json').write_text(json.dumps(gates,indent=2,allow_nan=False),encoding='utf-8')
print(pd.DataFrame(summaries)[['lineage','analysis','eligible_patient_n','median_C_total','median_C_between','median_C_within','lineage_gate']].to_string(index=False));print(json.dumps(gates,indent=2))
