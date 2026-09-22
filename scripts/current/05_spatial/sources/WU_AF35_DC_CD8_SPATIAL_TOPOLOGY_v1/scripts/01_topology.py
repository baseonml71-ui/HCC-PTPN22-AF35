from common import *
import numpy as np,itertools
from scipy.spatial import cKDTree

lock=verify_lock()
assert not (OUT/'QA/OUTCOME_START.json').exists(),'Do not silently rerun frozen topology analysis'
write_json('QA/OUTCOME_START.json',{'started_utc':now(),'freeze_time_utc':lock['freeze_time_utc'],'freeze_lock_sha256':sha(OUT/'WU_AF35_DC_CD8_TOPOLOGY_FREEZE_LOCK.json'),'protocol_sha256':lock['protocol_sha256']})
(OUT/'data/permutation_allocations').mkdir(exist_ok=True)
cells=read(OUT/'data/FROZEN_DC_CD8_HIGH_LOW_CELLS.tsv.gz')
elig=read(OUT/'WU_AF35_DC_CD8_TOPOLOGY_ELIGIBILITY.tsv')
elig['patient_order']=elig.patient_id.astype(int);elig['stratum_order']=elig.stratum.map({'LOW':0,'HIGH':1})
elig=elig.sort_values(['patient_order','region_id','stratum_order']).drop(columns=['patient_order','stratum_order']).reset_index(drop=True)
region_rows=[];null_rows=[];integrity=[];null_summary=[]
def nnd(xy,is_dc):return float(np.median(cKDTree(xy[~is_dc]).query(xy[is_dc],k=1,workers=1)[0]))
for gi,r in enumerate(elig.itertuples()):
    d=cells[(cells.region_id==r.region_id)&(cells.stratum==r.stratum)].sort_values('cell_id')
    xy=d[['x','y']].to_numpy(float);xy.flags.writeable=False
    observed_labels=(d.source_native_cell_type==DC).to_numpy()
    ndc=int(observed_labels.sum());ncd8=int((~observed_labels).sum())
    base=dict(group_index=gi,patient_id=r.patient_id,region_id=r.region_id,state='Pre',stratum=r.stratum,n_DC=ndc,n_CD8=ncd8,n_union=len(d))
    if ndc==0 or ncd8==0:
        region_rows.append(dict(**base,Observed_NND=np.nan,NullMedian_NND=np.nan,DC_CD8_PROXIMITY_ENRICHMENT=np.nan,status=r.count_reason,permutations=0));continue
    observed=nnd(xy,observed_labels)
    rng=np.random.Generator(np.random.PCG64(np.random.SeedSequence([lock['permutation_seed'],gi])))
    alloc=np.empty((lock['permutations'],ndc),dtype=np.int32)
    null=np.empty(lock['permutations']);coord_hash=hashlib.sha256(xy.tobytes()).hexdigest()
    for k in range(lock['permutations']):
        ids=rng.choice(len(d),size=ndc,replace=False);mask=np.zeros(len(d),bool);mask[ids]=True
        assert int(mask.sum())==ndc and int((~mask).sum())==ncd8
        assert hashlib.sha256(xy.tobytes()).hexdigest()==coord_hash
        alloc[k]=ids;null[k]=nnd(xy,mask)
        null_rows.append(dict(group_index=gi,permutation_index=k,patient_id=r.patient_id,region_id=r.region_id,stratum=r.stratum,permuted_NND=null[k],n_DC=ndc,n_CD8=ncd8,coordinates_sha256=coord_hash))
    path=OUT/f'data/permutation_allocations/group_{gi:02d}.npz'
    np.savez_compressed(path,dc_indices=alloc,cell_ids=d.cell_id.to_numpy(dtype='U'),xy=xy,observed_is_dc=observed_labels)
    nullmedian=float(np.median(null))
    valid=np.isfinite(observed) and np.isfinite(nullmedian) and observed>0 and nullmedian>0
    score=float(np.log2(nullmedian/observed)) if valid else np.nan
    region_rows.append(dict(**base,Observed_NND=observed,NullMedian_NND=nullmedian,DC_CD8_PROXIMITY_ENRICHMENT=score,status='EVALUABLE' if valid else 'NONPOSITIVE_OR_NONFINITE_NND',permutations=lock['permutations']))
    null_summary.append(dict(**base,Observed_NND=observed,NullMedian_NND=nullmedian,null_mean=float(null.mean()),null_sd=float(null.std(ddof=1)),null_q025=float(np.quantile(null,.025)),null_q975=float(np.quantile(null,.975)),DC_CD8_PROXIMITY_ENRICHMENT=score,permutations=lock['permutations'],seed_root=lock['permutation_seed'],group_stream=gi,coordinates_sha256=coord_hash,allocation_file=str(path.relative_to(OUT))))
    integrity.append(dict(group_index=gi,permutations=lock['permutations'],n_DC_preserved=True,n_CD8_preserved=True,coordinates_unchanged=True,xy_sha256=coord_hash,allocation_sha256=sha(path)))
    print(f'Completed frozen group {gi+1}/{len(elig)}; {lock["permutations"]} permutations',flush=True)
regions=pd.DataFrame(region_rows);save(regions,'WU_AF35_DC_CD8_TOPOLOGY_REGION_METRICS.tsv')
save(pd.DataFrame(null_rows),'data/WU_ALL_MARK_PERMUTATION_VALUES.tsv.gz')
save(pd.DataFrame(null_summary),'WU_AF35_DC_CD8_MARK_PERMUTATION_SUMMARY.tsv')
save(pd.DataFrame(integrity),'QA/PERMUTATION_INTEGRITY.tsv')
patient_rows=[];raw_rows=[]
for patient in lock['patients']:
    d=regions[(regions.patient_id==patient)&(regions.status=='EVALUABLE')]
    high=d[d.stratum=='HIGH'];low=d[d.stratum=='LOW']
    valid=len(high)>0 and len(low)>0
    h=high.DC_CD8_PROXIMITY_ENRICHMENT.mean();l=low.DC_CD8_PROXIMITY_ENRICHMENT.mean()
    hr=high.Observed_NND.mean();lr=low.Observed_NND.mean()
    patient_rows.append(dict(patient_id=patient,state='Pre',n_HIGH_regions=len(high),n_LOW_regions=len(low),HIGH_regions=';'.join(high.region_id),LOW_regions=';'.join(low.region_id),HIGH_proximity_enrichment=h,LOW_proximity_enrichment=l,DELTA_TOPOLOGY=h-l if valid else np.nan,status='EVALUABLE' if valid else 'MISSING_HIGH_AND_LOW' if not len(high) and not len(low) else 'MISSING_HIGH' if not len(high) else 'MISSING_LOW',weighting='equal_region_within_patient;equal_patient_between_patients'))
    raw_rows.append(dict(patient_id=patient,state='Pre',HIGH_DC_to_CD8_NND=hr,LOW_DC_to_CD8_NND=lr,HIGH_minus_LOW_raw_NND=hr-lr if valid else np.nan,HIGH_is_closer=bool(hr<lr) if valid else None,n_HIGH_regions=len(high),n_LOW_regions=len(low),distance_unit='source_CODEX_XY_units',role='SECONDARY_DESCRIPTIVE_ABUNDANCE_SENSITIVE'))
patients=pd.DataFrame(patient_rows);save(patients,'WU_AF35_DC_CD8_TOPOLOGY_PATIENT_METRICS.tsv');save(pd.DataFrame(raw_rows),'WU_AF35_DC_CD8_SECONDARY_RAW_NND.tsv')
v=patients[patients.status=='EVALUABLE'].DELTA_TOPOLOGY.to_numpy();N=len(v)
result=dict(N_evaluable=N,N_authoritative=len(lock['patients']),mean_DELTA_TOPOLOGY=float(v.mean()) if N else np.nan,median_DELTA_TOPOLOGY=float(np.median(v)) if N else np.nan,n_positive=int((v>0).sum()),positive_fraction=float((v>0).mean()) if N else np.nan,exact_two_sided_P=np.nan,bootstrap_CI_low=np.nan,bootstrap_CI_high=np.nan,bootstrap_resamples=0,sign_configurations=0)
if N<8:
    gate='NOT_EVALUABLE_LOW_PATIENT_COVERAGE'
else:
    # Bits enumerate all signs, including zeros without dropping patients.
    signs=2*((np.arange(2**N)[:,None]>>np.arange(N))&1)-1
    perm_mean=signs@v/N;obs=v.mean();p=float(np.mean(np.abs(perm_mean)>=abs(obs)-1e-12))
    save(pd.DataFrame({'configuration_index':np.arange(2**N),'signed_mean':perm_mean}),'data/WU_EXACT_PATIENT_SIGN_CONFIGURATIONS.tsv')
    rng=np.random.Generator(np.random.PCG64(lock['bootstrap_seed']))
    indices=rng.integers(0,N,size=(lock['bootstrap_n'],N));boot=v[indices].mean(axis=1)
    np.save(OUT/'data/BOOTSTRAP_PATIENT_INDICES.npy',indices)
    save(pd.DataFrame(indices+1), 'data/BOOTSTRAP_PATIENT_INDICES_1BASED.tsv')
    save(pd.DataFrame({'draw':np.arange(lock['bootstrap_n']),'mean_DELTA_TOPOLOGY':boot}),'data/WU_PATIENT_BOOTSTRAP_MEANS.tsv')
    ci=np.quantile(boot,[.025,.975],method='linear')
    positive=int((v>0).sum());consistent=positive>=int(np.ceil(.8*N))
    gate='PASS' if obs>0 and consistent and p<.05 else 'DIRECTIONALLY_SUPPORTIVE_NOT_SIGNIFICANT' if obs>0 and consistent else 'NOT_SUPPORTED'
    result.update(exact_two_sided_P=p,bootstrap_CI_low=ci[0],bootstrap_CI_high=ci[1],bootstrap_resamples=lock['bootstrap_n'],sign_configurations=2**N,positive_count_required=int(np.ceil(.8*N)),direction_consistency_pass=consistent,mean_positive=obs>0,primary_P_pass=p<.05,abundance_adjusted_positive=obs>0)
result['primary_gate']=gate;save(pd.DataFrame([result]),'WU_AF35_DC_CD8_PRIMARY_TEST.tsv')
save(pd.DataFrame([dict(endpoint='DC_CD8_ADJACENCY_SECONDARY',status='NOT_RUN_NO_FROZEN_GRAPH',reason='No frozen source-native CODEX cell adjacency graph in bounded source system; no graph created')]),'WU_AF35_DC_CD8_ADJACENCY_SENSITIVITY.tsv')
write_json('QA/OUTCOME_FINISH.json',{'finished_utc':now(),'primary_gate':gate,'frozen_lock_unchanged':sha(OUT/'WU_AF35_DC_CD8_TOPOLOGY_FREEZE_LOCK.json')==json.loads((OUT/'QA/OUTCOME_START.json').read_text())['freeze_lock_sha256']})
print(pd.DataFrame([result]).to_string(index=False),flush=True)
