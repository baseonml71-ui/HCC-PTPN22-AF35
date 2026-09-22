import importlib,sys
import numpy as np,pandas as pd
from scipy.stats import spearmanr
f=importlib.import_module('00_initialize');OUT=f.OUT
sys.path.insert(0,str(f.P/'10_software/phase3_py'));import h5py
lock=f.jl(OUT/'AF35_FUNCTIONAL_STATE_FREEZE_LOCK.json')
for n,h in lock['hashes'].items():assert f.sha(OUT/n)==h
genes=f.read(OUT/'inputs/AF35_ORDERED_GENES.tsv').gene.tolist();membership=f.read(OUT/'inputs/GENE_LAYER_MEMBERSHIP.tsv')
def cor(a,b):return float(np.corrcoef(a,b)[0,1])
def sp(a,b):return float(spearmanr(a,b).statistic)
network=[];coupling=[];layercoupling=[];mats={};curves={}
for c in f.COHORTS:
 d=f.read(f.H/f'inputs/{c}_GENE_Z_RESPONSE_BLIND.tsv',usecols=['patient_id']+genes).set_index('patient_id')
 ids=f.read(f.H/f'inputs/{c}_PATIENT_IDS.tsv').patient_id.tolist();assert d.index.tolist()==ids
 anchor=f.read(f.A/f'{c}_GENE_Z_AND_FULL_SCORES_RESPONSE_BLIND.tsv',usecols=['patient_id','PTPN22']).set_index('patient_id').loc[ids,'PTPN22']
 x=d[genes].to_numpy();r=d[genes].corr(method='spearman').to_numpy();mats[c]=r
 old=f.read(f.H/('AF35_PROGRAM_COHERENCE_GENE_CORRELATION_'+c+'.tsv'));assert np.max(abs(r-old[genes].to_numpy()))<1e-10
 cp=[]
 for i,g in enumerate(genes):
  edges=np.delete(r[i],i);network.append(dict(cohort=c,n=len(d),gene_rank=i+1,gene=g,gene_to_rest_rho=sp(x[:,i],(x.sum(axis=1)-x[:,i])/34),mean_connectivity=edges.mean(),median_connectivity=np.median(edges),positive_edge_fraction=(edges>0).mean()))
  cp.append(dict(cohort=c,n=len(d),gene_rank=i+1,gene=g,PTPN22_spearman=sp(x[:,i],anchor),PTPN22_pearson=cor(x[:,i],anchor)))
 curves[c]=pd.DataFrame(cp);coupling+=cp
 for layer,md in membership.groupby('layer_id',sort=False):
  cc=curves[c].set_index('gene').loc[md.gene]
  layercoupling.append(dict(cohort=c,layer_id=layer,gene_n=len(cc),mean_spearman=cc.PTPN22_spearman.mean(),median_spearman=cc.PTPN22_spearman.median(),mean_pearson=cc.PTPN22_pearson.mean(),positive_spearman_n=int((cc.PTPN22_spearman>0).sum())))
 f.save(d.reset_index().assign(PTPN22=anchor.to_numpy()),f'inputs/{c}_MOLECULAR_WITH_ANCHOR_RESPONSE_BLIND.tsv')
net=pd.DataFrame(network);cou=pd.DataFrame(coupling);cc=[]
for metric in ['mean_connectivity','median_connectivity','positive_edge_fraction','gene_to_rest_rho']:
 a=net[net.cohort==f.COHORTS[0]][metric];b=net[net.cohort==f.COHORTS[1]][metric]
 cc.append(dict(metric=metric,gene_n=35,Pearson=cor(a,b),Spearman=sp(np.round(a,12),np.round(b,12)),interpretation='DESCRIPTIVE_DEPENDENT_GENE_VECTOR'))
ccon=[]
for metric in ['PTPN22_spearman','PTPN22_pearson']:
 a=curves[f.COHORTS[0]][metric];b=curves[f.COHORTS[1]][metric]
 ccon.append(dict(coupling_vector=metric,gene_n=35,Pearson=cor(a,b),Spearman=sp(np.round(a,12),np.round(b,12))))
signs=cou.pivot(index='gene',columns='cohort',values='PTPN22_spearman');consistent=(np.sign(signs.iloc[:,0])==np.sign(signs.iloc[:,1]))
cou['direction']=np.where(cou.PTPN22_spearman>0,'POSITIVE',np.where(cou.PTPN22_spearman<0,'NEGATIVE','ZERO'));cou['cross_cohort_direction_consistent']=cou.gene.map(consistent)
for name,x in [('AF35_INTERNAL_NETWORK_GENE_METRICS.tsv',net),('AF35_INTERNAL_NETWORK_CROSS_COHORT.tsv',cc),('AF35_PTPN22_GENE_COUPLING.tsv',cou),('AF35_PTPN22_COUPLING_CONSERVATION.tsv',ccon),('AF35_PTPN22_FUNCTIONAL_LAYER_COUPLING.tsv',layercoupling)]:f.save(x,name)
print('I molecular network and anchor coupling complete',flush=True)
meta=f.read(OUT/'inputs/TCR_CELL_METADATA_RESPONSE_BLIND.tsv');states=f.read(OUT/'inputs/RELEASED_STATE_ORDER.tsv').released_state.tolist()
npz=np.load(f.C/'authority/03_processed_data/phase19a_af35_clonotype/GSE235863_AF35_RESPONSE_BLIND_PROJECTION.npz')
assert npz['barcode'].tolist()==meta.barcode.tolist() and npz['af35_gene_order'].tolist()==genes
path=f.P/'03_processed_data/phase3/GSE235863_nine_patients_scRNAseq_cd45_raw_counts.h5ad'
def decode(v):return np.array([s.decode() if isinstance(s,bytes) else str(s) for s in v])
with h5py.File(path,'r') as h:
 hg=decode(h['var/_index'][:]);barcodes=decode(h['obs/_index'][:]);total=np.asarray(h['obs/total_counts'][:],dtype=float)
 rows=meta.obs_row.to_numpy(dtype=int)-1
 # The mapping field is audited against barcodes, never assumed to use a given origin.
 if not np.array_equal(barcodes[rows],meta.barcode.to_numpy()):rows=meta.obs_row.to_numpy(dtype=int)
 assert np.array_equal(barcodes[rows],meta.barcode.to_numpy())
 gi=np.array([np.flatnonzero(hg==g).item() for g in genes]);x=h['X'];ptr=x['indptr'][:]
 rowmap=np.full(len(barcodes),-1,dtype=int);rowmap[rows]=np.arange(len(rows));gmap=np.full(len(hg),-1,dtype=int);gmap[gi]=np.arange(35)
 counts=np.zeros((len(rows),35),dtype=np.int64);nnz=len(x['indices'])
 for st in range(0,nnz,5000000):
  end=min(st+5000000,nnz);ix=x['indices'][st:end];gc=gmap[ix];keep=gc>=0
  pos=st+np.flatnonzero(keep);rr=rowmap[np.searchsorted(ptr[1:],pos,side='right')];valid=rr>=0
  raw=x['data'][st:end][keep][valid];assert np.all(raw==np.floor(raw))
  counts[rr[valid],gc[keep][valid]]=raw.astype(np.int64)
  if st==0 or end==nnz:print('H5AD CSR extraction',end,'/',nnz,flush=True)
 totals=total[rows];assert np.all(totals>0)
 # Independent direct row-slice extraction, different access pattern from global CSR scan.
 second=np.zeros_like(counts);totals_check=np.zeros(len(rows))
 for j,row in enumerate(rows):
  a,b=int(ptr[row]),int(ptr[row+1]);ix=x['indices'][a:b];vv=x['data'][a:b];gg=gmap[ix];keep=gg>=0
  second[j,gg[keep]]=vv[keep].astype(np.int64);totals_check[j]=vv.sum(dtype=np.float64)
 assert np.array_equal(second,counts) and np.array_equal(totals_check,totals)
f.save([dict(check='58872_by_35_raw_UMI_two_access_patterns',status='PASS',max_absolute_difference=int(abs(second-counts).max())),dict(check='all_cell_total_UMI_recomputed',status='PASS',max_absolute_difference=float(abs(totals_check-totals).max()))],'QA/H5AD_INDEPENDENT_EXTRACTION_VERIFICATION.tsv')
log=np.log1p(counts*(10000/totals[:,None]));score=((log-npz['gene_means'])/npz['gene_sds']).mean(axis=1)
err=float(abs(score-npz['af35_primary']).max());assert err<1e-7
f.save([dict(check='frozen_cell_AF35_reproduction',n=58872,max_absolute_error=err,status='PASS')],'QA/FROZEN_SINGLE_CELL_AF35_REPRODUCTION.tsv')
cells=pd.DataFrame(counts,columns=genes);cells.insert(0,'barcode',meta.barcode);cells['total_UMI']=totals.astype(np.int64);cells['frozen_AF35']=npz['af35_primary']
cells.to_csv(OUT/'inputs/CELL_AF35_RAW_COUNTS_AND_FROZEN_SCORE.tsv.gz',sep='\t',index=False)
records=[];score_records=[]
for scope in ['ALL_MAPPED','pre_P','pre_T','post_P','post_T']:
 mm=meta if scope=='ALL_MAPPED' else meta[(meta.timepoint+'_'+meta.tissue)==scope]
 for (patient,state),group in mm.groupby(['patient','sub_cluster'],sort=False):
  ind=group.index.to_numpy();cnt=counts[ind].sum(axis=0);den=totals[ind].sum();pb=np.log1p(10000*cnt/den)
  for j,g in enumerate(genes):records.append(dict(scope=scope,patient=patient,released_state=state,gene_rank=j+1,gene=g,cell_n=len(ind),total_UMI=int(den),gene_UMI=int(cnt[j]),pseudobulk_log1p_CP10K=pb[j],detected_cells=int((counts[ind,j]>0).sum()),cell_detection_fraction=float((counts[ind,j]>0).mean()),feature_present=True))
  score_records.append(dict(scope=scope,patient=patient,released_state=state,cell_n=len(ind),frozen_AF35_mean=float(npz['af35_primary'][ind].mean())))
ps=pd.DataFrame(records)
ps['within_patient_state_z']=ps.groupby(['scope','patient','gene']).pseudobulk_log1p_CP10K.transform(lambda v:(v-v.mean())/v.std(ddof=1) if v.std(ddof=1)>0 else np.nan)
summary=ps.groupby(['scope','released_state','gene_rank','gene'],sort=False).agg(patient_n=('patient','nunique'),mean_expression=('pseudobulk_log1p_CP10K','mean'),median_expression=('pseudobulk_log1p_CP10K','median'),mean_state_z=('within_patient_state_z','mean'),median_state_z=('within_patient_state_z','median'),standardized_patient_n=('within_patient_state_z','count'),positive_expression_patient_n=('gene_UMI',lambda v:int((v>0).sum())),mean_cell_detection_fraction=('cell_detection_fraction','mean'),total_cells=('cell_n','sum'),min_cells_per_patient=('cell_n','min')).reset_index()
eligible=meta.groupby(meta.timepoint+'_'+meta.tissue).patient.nunique().to_dict();eligible['ALL_MAPPED']=meta.patient.nunique()
summary['eligible_patients_in_scope']=summary.scope.map(eligible);summary['absent_patient_state_n']=summary.eligible_patients_in_scope-summary.patient_n;summary['gene_present_in_feature_matrix']=True
sc=pd.DataFrame(score_records);scs=sc.groupby(['scope','released_state'],sort=False).agg(patient_n=('patient','nunique'),mean_frozen_AF35=('frozen_AF35_mean','mean'),median_frozen_AF35=('frozen_AF35_mean','median'),cell_n=('cell_n','sum')).reset_index()
layer=ps.merge(membership[['gene','layer_id']],on='gene',how='left').groupby(['scope','patient','released_state','layer_id'],sort=False).agg(layer_mean_state_z=('within_patient_state_z','mean'),usable_gene_n=('within_patient_state_z','count')).reset_index()
ls=layer.groupby(['scope','released_state','layer_id'],sort=False).agg(patient_n=('layer_mean_state_z','count'),mean_layer_state_z=('layer_mean_state_z','mean'),median_layer_state_z=('layer_mean_state_z','median'),minimum_usable_genes=('usable_gene_n','min')).reset_index()
for name,x in [('AF35_STATE_LOCALIZATION_PATIENT_STATE.tsv',ps),('AF35_STATE_LOCALIZATION_SUMMARY.tsv',summary),('AF35_STATE_FROZEN_SCORE_PATIENT_STATE.tsv',sc),('AF35_STATE_FROZEN_SCORE_SUMMARY.tsv',scs),('AF35_FUNCTIONAL_LAYER_PATIENT_STATE.tsv',layer),('AF35_FUNCTIONAL_LAYER_STATE_SUMMARY.tsv',ls)]:f.save(x,name)
f.js('QA/MOLECULAR_STATE_COMPLETION_LOCK.json',dict(utc=f.now(),outcome_loaded=False,patients=9,cells=58872,released_states=20,raw_extraction_two_methods_agree=True,hashes={p.name:f.sha(p) for p in OUT.glob('AF35_*.tsv')}))
print('State pseudobulk complete;',len(ps),'gene-level patient/state/context rows;',len(sc),'patient/state/context units',flush=True)
