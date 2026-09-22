"""仅执行元数据、矩阵完整性和冻结成员覆盖审计；没有评分或效应分析。"""
from pathlib import Path
from collections import Counter
import re,json,hashlib,shutil,datetime,gzip
import pandas as pd
import numpy as np
import openpyxl

ROOT=Path(__file__).resolve().parents[1]
PROJECT=ROOT.parent
RAW=ROOT/'01_source/raw'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def tsv(rows,path):
    frame=rows if isinstance(rows,pd.DataFrame) else pd.DataFrame(rows)
    frame.to_csv(ROOT/path,sep='\t',index=False,encoding='utf-8',lineterminator='\n')
def note(path,text): (ROOT/path).write_text(text,encoding='utf-8')

# 权威冻结定义以原始字节复制并验证既有清单；复制代码不执行代码。
authority=PROJECT/'PTPN22_SIGNALING_FEEDBACK_MECHANISM_v1'
program_source=authority/'01_source/PTPN22_MECHANISM_PROGRAMS_FROZEN.tsv'
old_manifest=pd.read_csv(authority/'10_qa/DELIVERABLE_MANIFEST_SHA256.tsv',sep='\t')
expected=old_manifest.loc[old_manifest.path=='01_source/PTPN22_MECHANISM_PROGRAMS_FROZEN.tsv','sha256'].item()
assert sha(program_source)==expected
sources=[program_source,
 PROJECT/'HCC_ICI_project/04_scripts/clinical_translation/01_phase9_clinical_translation.R',
 PROJECT/'HCC_ICI_project/04_scripts/phase19/03_phase19b_response_blind_scoring.R']
local=[]
for p in sources:
    dest=ROOT/'01_source/local_authority'/p.name
    initial=sha(p);shutil.copyfile(p,dest);assert initial==sha(dest)==sha(p)
    local.append(dict(original_path=str(p),copied_path=dest.relative_to(ROOT).as_posix(),
        sha256=initial,old_manifest_expected_sha256=expected if p==program_source else 'CURRENT_SOURCE_HASH',
        copy_matches=True))
tsv(local,'07_qa/LOCAL_AUTHORITY_HASH_VERIFICATION.tsv')
programs=pd.read_csv(program_source,sep='\t',keep_default_na=False)
assert not programs.duplicated(['program','gene']).any()
af=programs.loc[programs.program=='PTPN22_activation_feedback_35','gene'].tolist()
requested='CD3D CD3E CD3G CD247 LCK FYN ZAP70 LAT LCP2 TRAT1 NFKB1 REL RELA NFKBIA TNFAIP3 STAT1 IRF1 IFIT1 IFIT2 IFIT3 ISG15 MX1 OAS1 CXCL9 CXCL10 IFNG TNF CCL5 NKG7 PRF1 GZMB PDCD1 LAG3 TIGIT TOX'.split()
assert af==requested and 'PTPN22' not in af

# 同次来源检索的原始网页工具提取文本，不从其他聊天重构供者资料。
supp=(RAW/'Chemin_2018_SuppMat_web_extraction.txt').read_text(encoding='utf-8')
donor_rows=re.findall(r'L\d+@P1: (\d+) (CC|TT) (\d+) (M|F)\s*(?:\n|$)',supp)
sup=pd.DataFrame(donor_rows,columns=['HD_no','source_paper_genotype','age','sex'])
assert len(sup)==26 and sup.HD_no.nunique()==26
sup['individual_id']='individual'+sup.HD_no
sup['age']=sup.age.astype(int);sup['sex']=sup.sex.map({'M':'Male','F':'Female'})
sup['source_location']='Chemin 2018 Supporting Information Table 1; PDF page 2; extracted lines L20-L45'
tsv(sup,'01_source/Chemin_2018_HD_TABLE1_TRANSCRIPTION.tsv')
meta=pd.read_csv(ROOT/'01_source/GEO_METADATA_PARSED.tsv',sep='\t',keep_default_na=False)
meta=meta.rename(columns={'individual':'individual_id','in vitro t cell activation':'condition_raw',
 'rs2476601 (ptpn22 1858) genotype':'genotype_raw','Sex':'sex',
 'hla-drb1 genotype':'HLA_DRB1','cytomegalovirus status':'CMV_status'})
meta['age']=meta.age.str.extract(r'^(\d+) years$')[0].astype(int)
activation='anti-CD3/CD28 beads (1bead/cell, Dynabeads) for 16 hours'
assert set(meta.condition_raw)=={'No',activation}
meta['condition_standardized']=meta.condition_raw.map({'No':'naive',activation:'activated'})
meta['processed_matrix_column']=meta.title
meta['pair_id']=meta.individual_id
meta['pair_complete']=False;meta['metadata_conflict']='NONE';meta['notes']=''
missing={'','--','?'}
def normalized(v): return 'MISSING' if str(v) in missing else str(v)
for individual,g in meta.groupby('individual_id',sort=False):
    complete=len(g)==2 and set(g.condition_standardized)=={'naive','activated'}
    conflicts=[c for c in ['genotype_raw','age','sex','HLA_DRB1','CMV_status'] if len({normalized(v) for v in g[c]})!=1]
    meta.loc[g.index,'pair_complete']=complete
    meta.loc[g.index,'metadata_conflict']=';'.join(conflicts) if conflicts else 'NONE'
    if g.HLA_DRB1.nunique()>1 and 'HLA_DRB1' not in conflicts:
        meta.loc[g.index,'notes']='HLA 在 naive 中字段缺失、activated 中为 --；两者均未知，不是已知等位基因冲突'
assert meta.GSM.nunique()==len(meta) and not meta.duplicated(['individual_id','condition_standardized']).any()
assert meta.pair_complete.all() and set(meta.metadata_conflict)=={'NONE'}
donors=meta.drop_duplicates('individual_id').merge(sup,on='individual_id',suffixes=('','_paper'),validate='one_to_one')
assert len(donors)==meta.individual_id.nunique()
assert (donors.age==donors.age_paper).all() and (donors.sex==donors.sex_paper).all()
# 用逐人供者表确立 GEO-原文关系；参考链只做独立一致性验证。
cross=donors.groupby('genotype_raw').source_paper_genotype.agg(lambda z:sorted(set(z))).to_dict()
assert cross=={'AA':['TT'],'GG':['CC']}
rs=json.loads((RAW/'rs2476601_RefSNP.json').read_text())
placements={p['seq_id']:p for p in rs['primary_snapshot_data']['placements_with_allele']}
assert rs['refsnp_id']=='2476601'
assert placements['NM_015967.7']['placement_annot']['is_aln_opposite_orientation']
assert 'NM_015967.7:c.1858C>T' in [a['hgvs'] for a in placements['NM_015967.7']['alleles']]
assert 'NP_057051.4:p.Trp620Arg' in [a['hgvs'] for a in placements['NP_057051.4']['alleles']]
tsv([dict(seq_id=k,opposite_orientation=v['placement_annot']['is_aln_opposite_orientation'],
 hgvs=';'.join(a['hgvs'] for a in v['alleles']),
 inserted_alleles=';'.join(a['allele']['spdi']['inserted_sequence'] for a in v['alleles']))
 for k,v in placements.items() if k in ['NC_000001.11','NM_015967.7','NM_015967.8','NP_057051.4']],
 '02_genotype/RS2476601_REFERENCE_PLACEMENTS.tsv')
mapping=[]
for _,d in donors.sort_values('HD_no',key=lambda s:s.astype(int)).iterrows():
    mapping.append(dict(individual_id=d.individual_id,GEO_genotype_raw=d.genotype_raw,
      source_paper_genotype=d.source_paper_genotype,
      allele_orientation='GEO 基因组正链 G/A；PTPN22 转录本反链 C/T；逐供者表格与 dbSNP 一致',
      protein_allotype_if_resolved='R620/R620' if d.genotype_raw=='GG' else 'W620/W620',
      risk_status_if_resolved='homozygous_non_risk' if d.genotype_raw=='GG' else 'homozygous_risk',
      evidence_source=f'GEO individual linkage; Supporting Table 1 HD {d.HD_no}; source abstract 1858T; rs2476601 RefSNP',
      resolution_status='RESOLVED'))
tsv(mapping,'02_genotype/GSE94859_GENOTYPE_MAPPING.tsv')
counts=[];balance=[]
for raw,paper,category in [('GG','CC','homozygous_non_risk'),('AG','CT','heterozygous'),('AA','TT','homozygous_risk')]:
    g=donors[donors.genotype_raw==raw]
    counts.append(dict(GEO_genotype_raw=raw,source_paper_genotype=paper,biological_category=category,
      donor_n=len(g),paired_donor_n=len(g),sample_n=int((meta.genotype_raw==raw).sum()),
      presence='PRESENT' if len(g) else 'ABSENT_NOT_MERGED',genotype_model='NOT_SELECTED'))
    if not len(g): continue
    hla=g.HLA_DRB1.map(normalized);cmv=g.CMV_status.map(normalized)
    allele=Counter(a for x in hla if x!='MISSING' for a in x.split(' / '))
    broad=Counter(a[:3] for x in hla if x!='MISSING' for a in x.split(' / '))
    balance.append(dict(GEO_genotype_raw=raw,source_paper_genotype=paper,donor_n=len(g),
      age_median=float(g.age.median()),age_min=int(g.age.min()),age_max=int(g.age.max()),age_missing_n=0,
      male_n=int((g.sex=='Male').sum()),female_n=int((g.sex=='Female').sum()),sex_missing_n=0,
      CMV_positive_n=int((cmv=='+').sum()),CMV_negative_n=int((cmv=='-').sum()),CMV_missing_n=int((cmv=='MISSING').sum()),
      HLA_known_donor_n=int((hla!='MISSING').sum()),HLA_missing_donor_n=int((hla=='MISSING').sum()),
      HLA_raw_allele_counts=json.dumps(dict(sorted(allele.items()))),HLA_broad_family_allele_counts=json.dumps(dict(sorted(broad.items()))),
      complete_separation='年龄范围重叠；性别及已知 CMV 两类均跨组；稀有 HLA 类别单组出现，HLA 不宜饱和调整',
      notes='描述性平衡审计；不做假设检验、结局模型、自动协变量选择或缺失值填补'))
tsv(counts,'02_genotype/GSE94859_GENOTYPE_COUNTS.tsv')
tsv(balance,'02_genotype/GSE94859_GENOTYPE_COVARIATE_BALANCE.tsv')

# 只读取工作簿。绝不保存回工作簿；校验读取前后原始字节未改变。
matrix_path=RAW/'GSE94859_RPKMs.xlsx';matrix_hash=sha(matrix_path)
book=openpyxl.load_workbook(matrix_path,read_only=True,data_only=True)
ws=book.active;header=next(ws.values);sheetname=ws.title
assert len(header)==len(set(header)) and len(book.sheetnames)==1
book.close()
matrix=pd.read_excel(matrix_path,engine='openpyxl')
x=matrix.iloc[:,2:].to_numpy(float)
assert list(matrix.columns[:2])==['Gene_symbol','RefSeq_ID']
assert set(matrix.columns[2:])==set(meta.processed_matrix_column)
assert len(matrix.columns[2:])==len(meta)
assert np.isfinite(x).all() and (x>=0).all()
ena=pd.read_csv(RAW/'SRP099568_ENA_run_metadata.tsv',sep='\t',keep_default_na=False)
assert set(ena.experiment_accession)==set(meta.SRA_accession)
linked=meta.merge(ena,left_on='SRA_accession',right_on='experiment_accession',validate='one_to_one')
assert (linked.biosample==linked.sample_accession).all()
tsv(linked[['GSM','individual_id','condition_standardized','genotype_raw','SRA_accession','biosample','run_accession',
 'library_name','library_layout','instrument_model','read_count','base_count','first_public','last_updated']],
 '01_source/GSE94859_SRA_SAMPLE_MAPPING.tsv')
required=['GSM','individual_id','condition_raw','condition_standardized','genotype_raw','age','sex','HLA_DRB1',
 'CMV_status','SRA_accession','processed_matrix_column','pair_id','pair_complete','metadata_conflict','notes']
tsv(meta[required],'01_source/GSE94859_SAMPLE_PAIRING_LEDGER.tsv')
# 第二份 GEO 元数据格式进行独立标题/基因型交叉检查。
series=gzip.open(RAW/'GSE94859_series_matrix.txt.gz','rt').read()
def series_values(key):
    line=next(z for z in series.splitlines() if z.startswith(key+'\t'))
    return [z.strip('"') for z in line.split('\t')[1:]]
assert series_values('!Sample_geo_accession')==meta.GSM.tolist()
assert series_values('!Sample_title')==meta.title.tolist()
genotype_line=next(z for z in series.splitlines() if z.startswith('!Sample_characteristics_ch1') and 'rs2476601' in z)
assert [z.strip('"').rsplit(': ',1)[-1] for z in genotype_line.split('\t')[1:]]==meta.genotype_raw.tolist()

symbol=matrix.Gene_symbol
dup=symbol.duplicated(keep=False);dates=symbol.map(lambda z:isinstance(z,datetime.datetime))
compound=symbol.map(lambda z:isinstance(z,str) and '|' in z)
anomaly=[]
for i in matrix.index[dup|dates|compound]:
    anomaly.append(dict(excel_row=int(i+2),gene_identifier_raw=str(symbol[i]),RefSeq_ID=matrix.RefSeq_ID[i],
      duplicate_identifier=bool(dup[i]),excel_date_identifier=bool(dates[i]),compound_identifier=bool(compound[i]),
      handling_phase0='RETAIN_UNCHANGED_NO_AGGREGATION',
      phase1_proposal='目标精确唯一行专用；非目标歧义行不进入评分；未经 RefSeq 溯源不得按日期反猜符号或拆分复合标识符'))
tsv(anomaly,'03_expression/GSE94859_IDENTIFIER_ANOMALIES.tsv')
sample_qc=[dict(processed_matrix_column=c,finite_n=int(np.isfinite(x[:,j]).sum()),missing_n=int(np.isnan(x[:,j]).sum()),
 negative_n=int((x[:,j]<0).sum()),nonzero_n=int((x[:,j]>0).sum()),GSM=meta.loc[meta.title==c,'GSM'].item()) for j,c in enumerate(matrix.columns[2:])]
tsv(sample_qc,'03_expression/GSE94859_SAMPLE_MATRIX_QC.tsv')
audit=[]
def metric(k,v,status='PASS',notes=''):audit.append(dict(check=k,value=v,status=status,notes=notes))
metric('expression_rows',len(matrix));metric('expression_sample_columns',x.shape[1]);metric('annotation_columns',2)
metric('worksheet',sheetname);metric('gene_identifier_system','Gene_symbol + pipe-delimited RefSeq_ID','PARTIAL','RefSeq 2014-12-21; hg38；包含复合符号和 Excel 日期型标识符')
metric('sample_identifier_system','GEO Sample_title',notes='32/32 标题、GSM、individual、SRX、BioSample 和 SRR 一一对应')
metric('author_expression_scale','FPKM in GEO processing; RPKM in filename/sheet','PARTIAL','作者连续非负长度/文库归一化数据；不是原始 count，不假称 TPM；命名差异保留')
metric('duplicate_sample_column_n',0);metric('missing_expression_values',int(np.isnan(x).sum()))
metric('negative_expression_values',int((x<0).sum()));metric('nonfinite_expression_values',int((~np.isfinite(x)).sum()))
metric('all_zero_rows',int((x==0).all(axis=1).sum()),notes='全矩阵 QC；没有按基因型或条件过滤')
metric('unique_identifier_values_including_dates',symbol.nunique(),'PARTIAL')
metric('duplicate_identifier_keys_including_dates',symbol[dup].nunique(),'PARTIAL','原始重复标识符，包含日期型值；不是全部已验证的不同生物学基因')
metric('duplicate_identifier_rows_including_dates',int(dup.sum()),'PARTIAL')
metric('duplicate_string_symbol_keys',symbol[dup&~dates].nunique(),'PARTIAL')
metric('excel_date_identifier_rows',int(dates.sum()),'PARTIAL','29 行原工作簿含日期型单元格；不涉及 PTPN22、AF35 或冻结程序精确目标行')
metric('compound_symbol_rows',int(compound.sum()),'PARTIAL','不拆分、不重复计入基因集')
metric('missing_gene_identifier_n',int(symbol.isna().sum()));metric('missing_RefSeq_ID_n',int(matrix.RefSeq_ID.isna().sum()))
metric('source_workbook_byte_unchanged',sha(matrix_path)==matrix_hash)
metric('duplicate_handling','NO_AGGREGATION','PARTIAL','不直接套用历史 bulk 重复符号均值；当前目标行都唯一，未来仅取冻结目标精确行')
metric('GSE94859_EXPRESSION_GATE','PARTIAL','PARTIAL','目标专用分析完整可用；全矩阵注释不能视为完全清洁')
tsv(audit,'03_expression/GSE94859_EXPRESSION_MATRIX_AUDIT.tsv')

def coverage(gene):
    indexes=matrix.index[symbol==gene].tolist()
    v=x[indexes,:] if indexes else np.empty((0,x.shape[1]))
    compound_matches=[int(i+2) for i,z in symbol.items() if isinstance(z,str) and '|' in z and gene in z.split('|')]
    valid=len(indexes)==1 and np.isfinite(v).all() and (v>=0).all() and np.std(v)>0
    return dict(gene=gene,exact_row_n=len(indexes),excel_row=';'.join(str(i+2) for i in indexes),
      RefSeq_ID=';'.join(str(matrix.RefSeq_ID[i]) for i in indexes),present=bool(indexes),
      duplicate_exact_symbol=len(indexes)>1,additional_compound_rows=';'.join(map(str,compound_matches)),
      finite_sample_n=int(np.isfinite(v).sum()),nonzero_sample_n=int((v>0).sum()),
      variable_across_all_samples=bool(len(indexes)==1 and np.std(v)>0),valid_unique_measurable=bool(valid),
      status='PASS' if valid else 'UNRESOLVED')
afcov=[coverage(g) for g in af];tsv(afcov,'03_expression/GSE94859_AF35_COVERAGE.tsv')
pt=coverage('PTPN22');pv=x[matrix.index[symbol=='PTPN22'],:]
pt.update(QC_min_all_samples=float(pv.min()),QC_max_all_samples=float(pv.max()),group_comparison='NOT_PERFORMED')
tsv([pt],'03_expression/GSE94859_PTPN22_COVERAGE.tsv')
prog_cov=[];prog_gene=[]
for name,g in programs.groupby('program',sort=False):
    if name=='PTPN22_activation_feedback_35':continue
    c=[coverage(z) for z in g.gene]
    prog_gene.extend([dict(program=name,**z) for z in c])
    prog_cov.append(dict(program=name,gene_n_frozen=len(c),gene_n_present=sum(z['present'] for z in c),
      gene_n_valid_unique=sum(z['valid_unique_measurable'] for z in c),missing_genes=';'.join(z['gene'] for z in c if not z['present']),
      duplicate_genes=';'.join(z['gene'] for z in c if z['duplicate_exact_symbol']),
      nonvalid_genes=';'.join(z['gene'] for z in c if not z['valid_unique_measurable']),
      coverage_fraction=sum(z['present'] for z in c)/len(c),status='PASS' if all(z['valid_unique_measurable'] for z in c) else 'PARTIAL',
      definition_source_sha256=expected))
tsv(prog_cov,'03_expression/GSE94859_FROZEN_PROGRAM_COVERAGE.tsv')
tsv(prog_gene,'03_expression/GSE94859_FROZEN_PROGRAM_GENE_COVERAGE.tsv')
assert all(z['valid_unique_measurable'] for z in afcov+prog_gene+[pt])
assert not any(z['additional_compound_rows'] for z in afcov+prog_gene+[pt])
summary=dict(paired_donor_n=len(donors),sample_n=len(meta),matrix_rows=len(matrix),
 genotype_counts=counts,covariate_balance=balance,AF35_valid_n=sum(z['valid_unique_measurable'] for z in afcov),
 frozen_programs=prog_cov,target_unique_gene_n=programs.gene.nunique(),
 duplicate_identifier_keys=int(symbol[dup].nunique()),duplicate_identifier_rows=int(dup.sum()),
 excel_date_rows=int(dates.sum()),allzero_rows=int((x==0).all(axis=1).sum()),matrix_sha256=matrix_hash,
 PTPN22=pt,source_publication_fulltext='UNAVAILABLE_HTTP403',source_supplement='WEB_TEXT_EXTRACT_ARCHIVED_PDF_BINARY_UNAVAILABLE')
note('07_qa/PHASE0_NUMERIC_SUMMARY.json',json.dumps(summary,ensure_ascii=False,indent=2))
print(json.dumps({k:summary[k] for k in ['paired_donor_n','sample_n','AF35_valid_n','target_unique_gene_n','duplicate_identifier_keys','duplicate_identifier_rows','excel_date_rows']},ensure_ascii=False))
print(pd.DataFrame(balance).to_string(index=False))
print('Phase 0 source / pairing / orientation / coverage assertions passed. No scores or biological effects calculated.')
