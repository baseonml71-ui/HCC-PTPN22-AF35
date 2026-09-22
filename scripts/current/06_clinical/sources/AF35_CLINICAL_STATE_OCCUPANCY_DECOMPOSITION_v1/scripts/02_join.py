from common import *
lock=jl(OUT/(P+'_FREEZE_LOCK.json'));assert not lock['outcome_loaded'] and not (OUT/(P+'_OUTCOME_JOIN_TIMESTAMP.txt')).exists()
for n,h in lock['hashes'].items():assert sha(OUT/n)==h,n
for n,h in lock['source_hashes'].items():assert sha(n)==h,n
stamp=now();write(P+'_OUTCOME_JOIN_TIMESTAMP.txt',stamp);js('QA/OUTCOME_OPEN_EVENT.json',dict(utc=stamp,freeze_utc=lock['utc'],freeze_lock_sha256=sha(OUT/(P+'_FREEZE_LOCK.json')),all_frozen_hashes_verified=True))
assert stamp>lock['utc']
source=D/'authority/03_processed_data/phase3/GSE235863_BARCODE_CLONOTYPE_STATE_PTPN22.tsv.gz'
labels=read(source,usecols=['patient','response_group']).drop_duplicates();assert labels.patient.is_unique
mapping={'Responder':1,'Nonresponder':0,'R':1,'NR':0,'responder':1,'nonresponder':0}
assert set(labels.response_group.dropna()).issubset(mapping),'Unknown source-native response labels; recover source coding without guessing'
labels['response']=labels.response_group.map(mapping)
save(labels,'inputs/RESPONSE_LABEL_MAP_AFTER_FREEZE.tsv')
pre=read(OUT/'inputs/PRETREATMENT_ELIGIBLE_CELLS_RESPONSE_BLIND.tsv');pc=read(OUT/(P+'_PATIENT_COMPONENTS_RESPONSE_BLIND.tsv'));design=read(OUT/'AF35_STATE_OCCUPANCY_MODEL_IDENTIFIABILITY.tsv')
clinical=pc.merge(labels,on='patient',validate='many_to_one',how='left');save(clinical,'inputs/PATIENT_COMPONENTS_AFTER_RESPONSE_JOIN.tsv')
feasible=[]
for compartment,g in pre.groupby('compartment',sort=True):
 pat=g[['patient']].drop_duplicates().merge(labels,on='patient',how='left',validate='one_to_one');n=int(pat.response.notna().sum());nr=int((pat.response==1).sum());nn=int((pat.response==0).sum())
 d=design[design.compartment==compartment]
 status='PASS' if n>=7 and min(nr,nn)>=2 else 'NOT_EVALUABLE_LOW_PRETREATMENT_COVERAGE'
 if status=='PASS' and (len(d)!=1 or not bool(d.iloc[0].identifiable)):status='NOT_EVALUABLE_DISCONNECTED_STATE_DESIGN'
 feasible.append(dict(compartment=compartment,source_native_timepoint='pre',patient_n=pat.patient.nunique(),response_evaluable_n=n,R_n=nr,NR_n=nn,missing_response_n=int(pat.response.isna().sum()),cell_n=len(g),status=status))
save(feasible,P+'_FEASIBILITY.tsv')
write(P+'_OUTCOME_FIREWALL_AUDIT.md',f'''# K 结局防火墙审计

全局组件冻结：{lock['utc']}；首次响应列解析事件：{stamp}。两时点严格先后，事件文件还记录原冻结JSON的SHA256。接入前逐项复核全部{len(lock['hashes'])}个冻结文件和{len(lock['source_hashes'])}个来源哈希。
冻结阶段只有明确无结局I元数据/细胞分数及Clone-State指定无结局列解析；post仅用于原始样本清单的排除记录及既有评分身份核对，不进入K模型或患者组件。源response_group列首次解析在本脚本事件落盘后。标签来自同一原Clone-State映射，逐患者唯一；未使用响应选择细胞/状态/区室/阈值/模型。
R/NR人数门在组件冻结后评估；所有通过区室运行，T在响应接入前已因6人不达标而停止。cap=3仅保留2.25%细胞，依结局前固定的过量丢失标准记NOT_RUN_LOW_CELL_COVERAGE，不调整cap。先前项目背景不是分析者历史绝对盲法，本次顺序与输入隔离由代码、文件及哈希审计支持。
''')
print(labels.to_string(index=False));print(pd.DataFrame(feasible).to_string(index=False))
