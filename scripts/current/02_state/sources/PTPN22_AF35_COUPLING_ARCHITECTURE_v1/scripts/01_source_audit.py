from pathlib import Path
import sys,json,hashlib,datetime,shutil
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent;P=BASE/'HCC_ICI_project'
sys.dont_write_bytecode=True;sys.path.insert(0,str(P/'10_software/phase3_py'))
import h5py
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(p,x):pd.DataFrame(x).to_csv(ROOT/p,sep='\t',index=False,na_rep='NA')
for d in ['01_source','02_manifest','03_primary','04_residual','05_state_architecture','06_patient_level','07_robustness','08_figures','09_qa','10_review','11_checkpoint','data','logs']:(ROOT/d).mkdir(exist_ok=True)
pass  # Non-computational private authorization copy excluded.
I=BASE/'AF35_FUNCTIONAL_ROLE_STATE_ARCHITECTURE_v1'
sources=[P/'03_processed_data/phase3/GSE235863_nine_patients_scRNAseq_cd45_raw_counts.h5ad',P/'03_processed_data/phase3/GSE235863_BARCODE_CLONOTYPE_STATE_PTPN22.tsv.gz',P/'03_processed_data/phase19a_af35_clonotype/GSE235863_AF35_RESPONSE_BLIND_PROJECTION.npz',P/'01_metadata/phase19a_af35_clonotype/PHASE19A_AF35_SINGLE_CELL_PROJECTION_CONTRACT.md',P/'01_metadata/phase19a_af35_clonotype/PHASE19A_AF35_GENE_PROJECTION_PARAMETERS.tsv',P/'04_scripts/phase19/01_phase19a_response_blind_projection.py',I/'inputs/AF35_ORDERED_GENES.tsv',I/'inputs/RELEASED_STATE_ORDER.tsv',I/'inputs/TCR_CELL_METADATA_RESPONSE_BLIND.tsv',I/'inputs/CELL_AF35_RAW_COUNTS_AND_FROZEN_SCORE.tsv.gz',I/'COMPLETE_SHA256_MANIFEST.tsv',I/'AF35_FUNCTIONAL_STATE_PROTOCOL_FREEZE_v1.md',I/'AF35_FUNCTIONAL_STATE_SOURCE_BINDINGS.tsv']
for rel in ['HCC_ICI_LATEST_PROJECT_CHECKPOINT.md','HCC_ICI_checkpoint_2026-08-31_v20_CONTROLLED_MANUSCRIPT_CONSTRUCTION_FREEZE.md','HCC_ICI_MAIN_FIGURE_SOURCE_LEDGER.tsv','HCC_ICI_v20R2_MAIN_FIGURE_EXACT_SOURCE_MAP.tsv','绘图视觉/HCC_ICI_FIVE_FIGURE_FINAL_ASSEMBLY_PHASE_A_V1/FIVE_FIGURE_FINAL_AUTHORITY_MANIFEST_V1.md','绘图视觉/HCC_ICI_FIGURE1_V2_6_RECONSTRUCTION_v1/source_copies/FROZEN_PANEL_LEDGER.tsv','HCC_ICI_BIOLOGICAL_RESOLUTION_BR1_v1/scripts/module_c_build_figure1_state_identity.R','HCC_ICI_BIOLOGICAL_RESOLUTION_BR1_v1/data/FIGURE1_STATE_PROGRAM_ASSOCIATION_PLOT_READY.tsv']:
    sources.append(BASE/rel)
for name in ['GSE149614_selected_gene_screen.rds','GSE326201_selected_gene_cell_screen.rds','GSE206325_compact_selected_gene_object.rds']:sources.append(P/'03_processed_data'/name)
for rel in ['04_scripts/01_GSE149614/01_candidate_screen_GSE149614.R','04_scripts/03_replication_scRNA/01_candidate_screen_GSE326201.R','04_scripts/02_GSE206325/01_extract_compact_GSE206325.R','04_scripts/config.R']:sources.append(P/rel)
expected={str(P/'03_processed_data/phase3/GSE235863_nine_patients_scRNAseq_cd45_raw_counts.h5ad'):'a1b5ce8a7f44cc50cba2458090dc079d7548f5567d1d9ee5237b791753c8d9d4',str(P/'03_processed_data/phase3/GSE235863_BARCODE_CLONOTYPE_STATE_PTPN22.tsv.gz'):'b4b611caab15865cfbbc6ebe4cf96ed3b4ef3f3d3c6c6910a68962df9dd617f3',str(P/'03_processed_data/phase19a_af35_clonotype/GSE235863_AF35_RESPONSE_BLIND_PROJECTION.npz'):'af65cf7a060604a267729bc50bfd6e5434f6f3370595421e33122fa936679176'}
im=pd.read_csv(I/'COMPLETE_SHA256_MANIFEST.tsv',sep='\t')
for r in im.itertuples():expected[str(I/r.path)]=r.sha256
records=[]
for p in dict.fromkeys(sources):
    h=sha(p);exp=expected.get(str(p));assert exp is None or h==exp,str(p)
    records.append(dict(path=str(p),bytes=p.stat().st_size,sha256=h,historical_sha256=exp,hash_status='PASS' if exp else 'CURRENT_SNAPSHOT'))
save('01_source/SOURCE_SHA256.tsv',records)
meta=pd.read_csv(I/'inputs/TCR_CELL_METADATA_RESPONSE_BLIND.tsv',sep='\t')
npz=np.load(P/'03_processed_data/phase19a_af35_clonotype/GSE235863_AF35_RESPONSE_BLIND_PROJECTION.npz')
assert np.array_equal(meta.barcode,npz['barcode']) or set(meta.barcode)==set(npz['barcode'])
af=pd.read_csv(I/'inputs/AF35_ORDERED_GENES.tsv',sep='\t');assert af.gene.tolist()==npz['af35_gene_order'].tolist() and len(af)==35 and 'PTPN22' not in set(af.gene)
save('09_qa/AF35_MEMBERSHIP_AUDIT.tsv',af.assign(in_frozen_projection=True,usable=npz['gene_usable']))
with h5py.File(sources[0],'r') as h:
    genes=h['var/_index'].asstr()[:];assert all(g in genes for g in af.gene)
    assert 'PTPN22' in genes
    rows=meta.obs_row.to_numpy()-1;bar=h['obs/_index'].asstr()[:]
    assert np.array_equal(bar[rows],meta.barcode)
    dictionary=[]
    for name in ['major_cluster','sub_cluster','patient','sample','tissue']:
        x=h['obs'][name]
        if isinstance(x,h5py.Group):vals=x['categories'].asstr()[:];codes=x['codes'][:];a=vals[codes]
        else:a=x.asstr()[:]
        dictionary.append(dict(field=name,values=';'.join(map(str,np.unique(a))),n_unique=len(np.unique(a))))
    save('01_source/GSE235863_RELEASED_METADATA_DICTIONARY.tsv',dictionary)
    covs=pd.DataFrame(dict(barcode=meta.barcode,total_UMI=h['obs/total_counts'][:][rows],detected_genes=h['obs/n_genes_by_counts'][:][rows]))
    covs.to_csv(ROOT/'data/TECHNICAL_COVARIATES.tsv.gz',sep='\t',index=False)
cells=meta.groupby(['timepoint','tissue']).agg(patient_n=('patient','nunique'),cell_n=('barcode','size'),state_n=('sub_cluster','nunique')).reset_index();save('09_qa/CELL_COUNT_AUDIT.tsv',cells)
meta.groupby(['patient','timepoint','tissue','sub_cluster']).size().rename('cell_n').reset_index().to_csv(ROOT/'09_qa/PATIENT_STATE_COVERAGE_BEFORE_EFFECTS.tsv',sep='\t',index=False)
for name in ['AF35_ORDERED_GENES.tsv','RELEASED_STATE_ORDER.tsv','TCR_CELL_METADATA_RESPONSE_BLIND.tsv']:shutil.copyfile(I/'inputs'/name,ROOT/'01_source'/name)
audit=[]
rstruct=pd.read_csv(ROOT/'01_source/R_OBJECT_STRUCTURE_AUDIT.tsv',sep='\t')
for r in rstruct.itertuples():
    state='released res.3 cluster IDs; six functional-state modules are continuous, not disjoint labels' if r.source=='GSE149614' else 'broad compartment only; six continuous modules, no frozen disjoint T-state labels' if r.source=='GSE326201' else 'six Figure 1 released-annotation-derived states; original treated tumor universe'
    audit.append(dict(source=r.source,patient_n=r.patient_n,cell_n=r.cell_n,state_system=state,PTPN22_available=True,AF35_35of35=False,patient_ID_available=True,treatment_timepoint_information='source metadata; not used for outcome selection',normalization_authority='source normalized matrix' if r.source=='GSE149614' else 'frozen log1pCP10K from full library size',eligible_primary=False,eligible_secondary=False,reason=f'No existing source-specific frozen AF35 single-cell scoring; compact object retains {r.AF35_genes_present}/35 gene expression columns. No new scoring reconstruction authorized.',count_scope='all cells in compact source object, not eligible analysis T cells'))
audit.append(dict(source='GSE235863',patient_n=9,cell_n=58872,state_system='20 frozen TCR released sub_cluster labels; aliases not harmonized with H5AD',PTPN22_available=True,AF35_35of35=True,patient_ID_available=True,treatment_timepoint_information='pre/post x tumor(T)/blood(P); baseline tumor 6 patients and 9088 cells',normalization_authority='PTPN22 log1pCP10K; AF35 Phase19A fixed global gene means/sample SD and equal weight',eligible_primary='PENDING_USER_SOURCE_ADJUSTMENT',eligible_secondary=True,reason='Eligible exact frozen TCR-mapped source; not an independent Figure 1 source. Primary-source adjustment requested before effect calculation.',count_scope='all existing TCR-mapped universe; contexts never pooled for inference'))
save('01_source/PTPN22_AF35_COUPLING_SOURCE_AUDIT.tsv',audit)
protected=[]
for branch in ['PTPN22_STATE_GATING_MECHANISM_v1','PTPN22_SIGNALING_FEEDBACK_MECHANISM_v1']:
    for p in (BASE/branch).rglob('*'):
        if p.is_file():protected.append(dict(path=str(p),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns))
save('09_qa/CLOSED_BRANCH_FILE_BASELINE.tsv',protected)
print('SOURCE HASH AND IDENTITY AUDIT PASS',len(records),flush=True)
print(cells.to_string(index=False),flush=True)
