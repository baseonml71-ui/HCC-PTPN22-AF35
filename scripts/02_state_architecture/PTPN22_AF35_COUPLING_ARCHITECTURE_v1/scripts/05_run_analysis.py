from pathlib import Path
import sys, hashlib, json, datetime
import numpy as np
import pandas as pd
from coupling_math import decompose,correlation,technical_residual

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent;P=BASE/'HCC_ICI_project';SOURCE='GSE235863'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(path,rows):pd.DataFrame(rows).to_csv(ROOT/path,sep='\t',index=False,na_rep='NA',float_format='%.15g')
lock=json.loads((ROOT/'02_manifest/FREEZE_LOCK.json').read_text())
for path,key in [('02_manifest/PTPN22_AF35_COUPLING_ANALYSIS_MANIFEST_V1.md','manifest_sha256'),('01_source/PTPN22_AF35_COUPLING_SOURCE_AUDIT.tsv','source_audit_sha256'),('01_source/SOURCE_SHA256.tsv','source_hash_manifest_sha256')]:assert sha(ROOT/path)==lock[key]
assert not (ROOT/'03_primary/PTPN22_AF35_PATIENT_COUPLING_DECOMPOSITION.tsv').exists(),'Do not silently overwrite effects'
(ROOT/'logs/EFFECT_RUN_STARTED.json').write_text(json.dumps(dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),freeze_sha256=lock['manifest_sha256']),indent=2))
meta=pd.read_csv(ROOT/'01_source/TCR_CELL_METADATA_RESPONSE_BLIND.tsv',sep='\t')
states=pd.read_csv(ROOT/'01_source/RELEASED_STATE_ORDER.tsv',sep='\t').released_state.tolist()
npz=np.load(P/'03_processed_data/phase19a_af35_clonotype/GSE235863_AF35_RESPONSE_BLIND_PROJECTION.npz')
scores=pd.DataFrame(dict(barcode=npz['barcode'],AF35=npz['af35_primary']))
cols=['barcode','ptpn22_raw_count','ptpn22_log1p_cp10k','total_raw_umi']
expr=pd.read_csv(P/'03_processed_data/phase3/GSE235863_BARCODE_CLONOTYPE_STATE_PTPN22.tsv.gz',sep='\t',usecols=cols)
cov=pd.read_csv(ROOT/'data/TECHNICAL_COVARIATES.tsv.gz',sep='\t')
data=meta.merge(scores,on='barcode',validate='one_to_one').merge(expr,on='barcode',validate='one_to_one').merge(cov,on='barcode',validate='one_to_one')
assert len(data)==58872 and data.barcode.is_unique
assert np.allclose(data.total_raw_umi,data.total_UMI)
data=data.loc[(data.timepoint=='pre')&(data.tissue=='T')].copy()
assert len(data)==9088 and data.patient.nunique()==6 and set(data.sub_cluster)==set(states)
assert np.isfinite(data[['AF35','ptpn22_raw_count','ptpn22_log1p_cp10k','total_UMI','detected_genes']]).all().all()
assert np.allclose(data.ptpn22_log1p_cp10k,np.log1p(10000*data.ptpn22_raw_count/data.total_raw_umi))
save('data/ANALYSIS_INPUT_RESPONSE_BLIND.tsv.gz',data)
zero_fraction=float((data.ptpn22_raw_count==0).mean());detection_trigger=zero_fraction>=.50

def sign(x):return int(np.sign(x)) if np.isfinite(x) and abs(x)>1e-12 else 0
def summarize_one(frame,variant,x=None,y=None):
    x=frame.ptpn22_log1p_cp10k.to_numpy() if x is None else x
    y=frame.AF35.to_numpy() if y is None else y
    out=decompose(x,y,frame.sub_cluster.to_numpy())
    row=dict(source=SOURCE,patient=frame.patient.iloc[0] if len(frame) else '',analysis=variant,cell_n=len(frame),state_n=frame.sub_cluster.nunique(),PTPN22_variance=float(np.var(x)) if len(x) else np.nan,AF35_variance=float(np.var(y)) if len(y) else np.nan)
    if out is None:
        row.update({k:np.nan for k in ['C_total','C_between','C_within','identity_error','residual_Pearson','residual_Spearman','between_abs_gt_within','sign_total','sign_between','sign_within','between_share','within_share']});row.update(eligibility='NA_ZERO_VARIANCE_OR_INSUFFICIENT',geometry='NOT_EVALUABLE')
        return row,[]
    detail=out.pop('states');row.update(out);t,b,w=(out[k] for k in ['C_total','C_between','C_within'])
    geometry='TOTAL_NEAR_ZERO' if abs(t)<.02 else 'OPPOSING_COMPONENTS' if b*w<0 else 'CONCORDANT_POSITIVE' if t>0 else 'CONCORDANT_NEGATIVE'
    shares=abs(t)>=.02 and sign(b)==sign(t)==sign(w) and sign(t)!=0
    row.update(eligibility='ELIGIBLE',between_abs_gt_within=abs(b)>abs(w),sign_total=sign(t),sign_between=sign(b),sign_within=sign(w),geometry=geometry,between_share=b/t if shares else np.nan,within_share=w/t if shares else np.nan)
    return row,detail

primary=[];residual=[];centroids=[];specific=[];pseudobulk=[];rare=[];technical=[];identity=[]
for patient,f in data.groupby('patient',sort=True):
    row,detail=summarize_one(f,'PRIMARY');primary.append(row)
    residual.append({k:row[k] for k in ['source','patient','cell_n','state_n','eligibility','residual_Pearson','residual_Spearman']})
    lookup={r['state']:r for r in detail}
    pb=dict(source=SOURCE,patient=patient,cell_n=len(f),mean_PTPN22=f.ptpn22_log1p_cp10k.mean(),mean_AF35=f.AF35.mean())
    for state in states:
        g=f[f.sub_cluster==state];n=len(g);pb['occupancy_'+state]=n/len(f)
        cr=dict(source=SOURCE,patient=patient,state=state,cell_n=n,mean_PTPN22=g.ptpn22_log1p_cp10k.mean(),mean_AF35=g.AF35.mean())
        cr.update({k:lookup.get(state,{}).get(k,np.nan) for k in ['weight','mean_Xz','mean_Yz','Cov_state','between_contribution','within_contribution']});centroids.append(cr)
        eligible=n>=20 and g.ptpn22_log1p_cp10k.var(ddof=0)>0 and g.AF35.var(ddof=0)>0
        specific.append(dict(source=SOURCE,patient=patient,state=state,cell_n=n,eligible=eligible,Pearson=correlation(g.ptpn22_log1p_cp10k,g.AF35) if eligible else np.nan,Spearman=correlation(g.ptpn22_log1p_cp10k,g.AF35,True) if eligible else np.nan,reason='ELIGIBLE' if eligible else 'LT20' if n<20 else 'ZERO_VARIANCE'))
    pseudobulk.append(pb)
    keep=f.sub_cluster.map(f.sub_cluster.value_counts())>=20
    rr,_=summarize_one(f.loc[keep],'RARE_STATE_GE20');rr['patient']=patient;rr.update(original_cell_n=len(f),retained_fraction=float(keep.mean()));rare.append(rr)
    x,y,rank=technical_residual(f.ptpn22_log1p_cp10k,f.AF35,f.total_UMI,f.detected_genes)
    tr,_=summarize_one(f,'TECHNICAL_ADJUSTED',x,y);tr.update(design_rank=rank,PTPN22_zero_fraction=float((f.ptpn22_raw_count==0).mean()),detection_triggered=detection_trigger);technical.append(tr)
    if detection_trigger:
        dr,_=summarize_one(f,'DETECTION_BINARY',(f.ptpn22_raw_count.to_numpy()>0).astype(float));dr.update(design_rank=np.nan,PTPN22_zero_fraction=float((f.ptpn22_raw_count==0).mean()),detection_triggered=True);technical.append(dr)
    else:
        technical.append(dict(source=SOURCE,patient=patient,analysis='DETECTION_BINARY',eligibility='NOT_TRIGGERED',PTPN22_zero_fraction=float((f.ptpn22_raw_count==0).mean()),detection_triggered=False))

U='COUPLING_ARCHITECTURE_UNRESOLVED'
def core_gate(f,min_n=3):
    f=f[f.eligibility=='ELIGIBLE'];n=len(f)
    if n<min_n:return U
    t=f.C_total.median();sg=sign(t)
    if abs(t)<.02 or (f.C_total.map(sign)==sg).mean()<2/3:return U
    b=f.C_between.abs().median();w=f.C_within.abs().median()
    if b>=.02 and w>=.02 and (f.C_between*f.C_within<0).mean()>=2/3:return 'OPPOSING_COUPLING_ARCHITECTURE'
    if (f.C_between.abs()>f.C_within.abs()).mean()>=2/3 and b>=1.5*w and (f.residual_Pearson.abs().median()<.10 or w<=b/1.5):return 'STATE_SEPARATION_PREDOMINANT'
    rp=(f.residual_Pearson.map(sign)==sg).mean();rs=(f.residual_Spearman.map(sign)==sg).mean()
    if (f.C_within.abs()>f.C_between.abs()).mean()>=2/3 and w>=1.5*b and rp>=2/3 and rs>=2/3:return 'WITHIN_STATE_PREDOMINANT'
    if b>=.02 and w>=.02:return 'MIXED_COUPLING_ARCHITECTURE'
    return U

lopo=[];gate_rows=[]
def full_gate(rows,name):
    f=pd.DataFrame(rows);f=f[f.eligibility=='ELIGIBLE'];core=core_gate(f,4);med=f[['C_total','C_between','C_within']].median();stable=True;changes={c:[] for c in med.index}
    for patient in f.patient:
        g=f[f.patient!=patient];lm=g[med.index].median();lg=core_gate(g);gs=lg==core and sign(lm.C_total)==sign(med.C_total);stable=stable and gs
        for c in med.index:
            delta=float(lm[c]-med[c]);changes[c].append(abs(delta))
            lopo.append(dict(source=SOURCE,analysis=name,left_out_patient=patient,eligible_patient_n=len(g),metric=c,full_estimate=med[c],lopo_estimate=lm[c],change=delta,sign_stable=sign(lm[c])==sign(med[c]),full_core_gate=core,lopo_core_gate=lg,qualitative_stable=gs))
    final=core if stable else U
    d=dict(source=SOURCE,analysis=name,eligible_patient_n=len(f),core_gate=core,final_gate=final,LOPO_qualitative_stable=stable,median_C_total=med.C_total,median_C_between=med.C_between,median_C_within=med.C_within,median_abs_between=f.C_between.abs().median(),median_abs_within=f.C_within.abs().median(),between_dominant_fraction=float((f.C_between.abs()>f.C_within.abs()).mean()),within_dominant_fraction=float((f.C_within.abs()>f.C_between.abs()).mean()),total_direction_fraction=float((f.C_total.map(sign)==sign(med.C_total)).mean()),opposing_fraction=float((f.C_between*f.C_within<0).mean()),median_residual_Pearson=f.residual_Pearson.median(),median_residual_Spearman=f.residual_Spearman.median())
    for c,v in changes.items():d['LOPO_max_abs_change_'+c]=max(v) if v else np.nan
    gate_rows.append(d);return d

gmain=full_gate(primary,'PRIMARY');grare=full_gate(rare,'RARE_STATE_GE20');gtech=full_gate([r for r in technical if r['analysis']=='TECHNICAL_ADJUSTED'],'TECHNICAL_ADJUSTED')
if detection_trigger:full_gate([r for r in technical if r['analysis']=='DETECTION_BINARY'],'DETECTION_BINARY')
for r in primary+rare+technical:
    if r.get('eligibility')=='ELIGIBLE':identity.append({k:r[k] for k in ['source','patient','analysis','cell_n','C_total','C_between','C_within','identity_error']}|dict(identity_pass=abs(r['identity_error'])<1e-12))
save('03_primary/PTPN22_AF35_PATIENT_COUPLING_DECOMPOSITION.tsv',primary)
save('04_residual/PTPN22_AF35_STATE_RESIDUAL_COUPLING.tsv',residual)
save('05_state_architecture/PTPN22_AF35_STATE_CENTROIDS.tsv',centroids)
save('05_state_architecture/PTPN22_AF35_STATE_SPECIFIC_COUPLING.tsv',specific)
save('06_patient_level/PTPN22_AF35_PATIENT_PSEUDOBULK.tsv',pseudobulk)
save('07_robustness/PTPN22_AF35_RARE_STATE_SENSITIVITY.tsv',rare)
save('07_robustness/PTPN22_AF35_TECHNICAL_SENSITIVITY.tsv',technical)
save('07_robustness/PTPN22_AF35_COUPLING_LOPO.tsv',lopo)
save('09_qa/COUPLING_IDENTITY_AUDIT.tsv',identity)
save('07_robustness/SOURCE_GATE_SUMMARY.tsv',gate_rows)
cent=pd.DataFrame(centroids);sp=pd.DataFrame(specific);summ=[];ss=[]
for state in states:
    c=cent[(cent.state==state)&(cent.cell_n>0)];g=sp[(sp.state==state)&sp.eligible]
    summ.append(dict(source=SOURCE,state=state,patient_n=len(c),cell_n=int(c.cell_n.sum()),patient_equal_mean_PTPN22=c.mean_PTPN22.mean(),patient_equal_mean_AF35=c.mean_AF35.mean(),median_PTPN22=c.mean_PTPN22.median(),median_AF35=c.mean_AF35.median()))
    ss.append(dict(source=SOURCE,state=state,eligible_patient_n=len(g),median_Pearson=g.Pearson.median(),min_Pearson=g.Pearson.min(),max_Pearson=g.Pearson.max(),positive_n=int((g.Pearson>1e-12).sum()),negative_n=int((g.Pearson< -1e-12).sum()),zero_n=int((g.Pearson.abs()<=1e-12).sum()),median_Spearman=g.Spearman.median()))
save('05_state_architecture/STATE_CENTROID_PATIENT_EQUAL_SUMMARY.tsv',summ)
save('05_state_architecture/STATE_SPECIFIC_PATIENT_SUMMARY.tsv',ss)
pb=pd.DataFrame(pseudobulk);save('06_patient_level/PATIENT_OVERALL_DESCRIPTIVE_ASSOCIATION.tsv',[dict(source=SOURCE,patient_n=len(pb),Pearson=correlation(pb.mean_PTPN22,pb.mean_AF35),Spearman=correlation(pb.mean_PTPN22,pb.mean_AF35,True),scope='DESCRIPTIVE_ONLY_NO_OUTCOME')])
techmaterial=gtech['final_gate']!=gmain['final_gate'] or sign(gtech['median_C_total'])!=sign(gmain['median_C_total']) or gtech['final_gate']==U
rarechange=grare['final_gate']!=gmain['final_gate'] or sign(grare['median_C_total'])!=sign(gmain['median_C_total'])
increment='MODERATE' if gmain['final_gate']!=U and not techmaterial and not rarechange else 'NONE'
final=dict(primary_source=SOURCE,context='pre/T',patient_n=6,cell_n=9088,state_n=20,COUPLING_ARCHITECTURE=gmain['final_gate'],CROSS_SOURCE_REPLICATION='NOT_INTERPRETABLE',cross_source_reason='No independent existing source satisfies frozen source-specific AF35 single-cell scoring requirement',COUPLING_MANUSCRIPT_INCREMENT=increment,technical_sensitivity='TECHNICAL_SENSITIVITY_MATERIAL' if techmaterial else 'TECHNICAL_ARCHITECTURE_RETAINED',rare_state_architecture_changed=rarechange,PTPN22_zero_fraction=zero_fraction,detection_triggered=detection_trigger,primary_summary=gmain,technical_summary=gtech,rare_state_summary=grare,independent_replication_available=False,outcome_used=False,causal_claim=False,analysis_completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(ROOT/'11_checkpoint/FINAL_GATES.json').write_text(json.dumps(final,indent=2,allow_nan=False),encoding='utf-8')
print(json.dumps(final,indent=2))
