from common import *
import numpy as np,shutil
for name in ['authority','QA','data','figures','scripts']:(OUT/name).mkdir(exist_ok=True)
assert not (OUT/'WU_AF35_DC_CD8_TOPOLOGY_FREEZE_LOCK.json').exists(),'Do not overwrite an existing protocol freeze'
pass  # Non-computational private authorization copy excluded.
pass  # Non-computational private authorization copy excluded.
manifest=read(C/'COMPLETE_SHA256_MANIFEST.tsv')
ver=[]
for r in manifest.itertuples():
    p=C/r.path;h=sha(p);ver.append(dict(path=str(p),expected=r.sha256,observed=h,status='PASS' if h==r.sha256 else 'FAIL'))
save(pd.DataFrame(ver),'QA/ANALYSIS_C_MANIFEST_VERIFICATION.tsv')
assert all(r['status']=='PASS' for r in ver),'Analysis C authority changed'
source_roles={
 'eligible_manifests/WU_ELIGIBLE_MAPPED_CELL_MEMBERSHIP.tsv.gz':'cell XY, labels, IDs, patient, region, frozen stratum',
 'eligible_manifests/WU_ELIGIBLE_REGION_SPOTS.tsv':'frozen region-spot eligibility and residual mapping',
 'eligible_manifests/WU_ELIGIBLE_REGIONS.tsv':'17-region patient authority',
 'eligible_manifests/WU_ELIGIBLE_PATIENTS.tsv':'10-patient authority',
 'eligible_manifests/WU_FROZEN_SECTION_SPOTS.tsv':'frozen whole-section residual and HIGH/LOW strata',
 'eligible_manifests/WU_SECTION_QUARTILE_THRESHOLDS.tsv':'frozen type-7 thresholds',
 'eligible_manifests/WU_NATIVE_CELLTYPE_CODEBOOK.tsv':'source-native labels',
 'WU_AF35_CONTEXT_SOURCE_BINDING_AUDIT.md':'source-binding semantics',
 'WU_AF35_ANALYSIS_CONTRACT.md':'Analysis C counting and model rules',
 'COMPLETE_SHA256_MANIFEST.tsv':'Analysis C immutable delivery manifest'}
binding=[dict(path=str(C/p),role=role,bytes=(C/p).stat().st_size,sha256=sha(C/p),hashed_at_utc=now()) for p,role in source_roles.items()]
save(pd.DataFrame(binding),'WU_AF35_DC_CD8_TOPOLOGY_SOURCE_BINDING_AUDIT.tsv')
# Protection extends beyond the inputs; only this new module is writable.
scope_dirs=[C,ROOT/'绘图视觉',ROOT/'HCC_ICI_project/06_figures',ROOT/'HCC_ICI_project/07_manuscript',ROOT/'写稿蓝图',ROOT/'_v20R2_source_completion']
scope_dirs += [ROOT/n for n in ['AF35_COMPONENT_ABLATION_ASSOCIATION_ROBUSTNESS_v1','AF35_MATCHED_RANDOM_CLINICAL_NULL_v1','CLINICAL_ANCHOR_READOUT_DECOMPOSITION_v1','GSE235863_CLONE_STATE_DECOMPOSITION_v1']]
protected=set()
for folder in scope_dirs:
    if folder.exists():protected.update(p for p in folder.rglob('*') if p.is_file())
old=json.loads((C/'QA/PROTECTED_INPUT_BASELINE.json').read_text(encoding='utf-8'))
protected.update(Path(r['path']) for r in old)
protected.update(p for p in ROOT.iterdir() if p.is_file() and ('manuscript' in p.name.lower() or 'legend' in p.name.lower()))
protected.add(ROOT/'HCC_ICI_project/01_metadata/clinical_translation/PTPN22_CLINICAL_AXES_FROZEN.tsv')
records=[dict(path=str(p),bytes=p.stat().st_size,sha256_before=sha(p)) for p in sorted(protected)]
write_json('QA/PROTECTED_BASELINE.json',{'created_utc':now(),'scope_directories':[str(p) for p in scope_dirs],'files':records})
print('Protected baseline',len(records),flush=True)
edges=read(C/'eligible_manifests/WU_ELIGIBLE_MAPPED_CELL_MEMBERSHIP.tsv.gz')
ctx=read(C/'eligible_manifests/WU_ELIGIBLE_REGION_SPOTS.tsv')
regions=read(C/'eligible_manifests/WU_ELIGIBLE_REGIONS.tsv')
patients=sorted(ctx.patient_id.unique(),key=int)
assert len(ctx)==5205 and len(regions)==17 and len(patients)==10 and len(edges)==151302
assert set(edges.state)=={'Pre'} and not edges.duplicated(['region_id','cell_id']).any()
assert np.isfinite(edges[['x','y']]).all().all()
join=edges.merge(ctx[['region_id','spot_id','patient_id','state','stratum']],on=['region_id','spot_id'],suffixes=('','_frozen'),validate='many_to_one')
assert len(join)==len(edges)
for col in ['patient_id','state','stratum']:assert np.array_equal(join[col],join[col+'_frozen'])
assert edges.groupby('region_id').patient_id.nunique().max()==1
pair=edges[edges.source_native_cell_type.isin([DC,CD8]) & edges.stratum.isin(['HIGH','LOW'])].copy()
pair=pair.sort_values(['patient_id','region_id','stratum','cell_id']).reset_index(drop=True)
save(pair,'data/FROZEN_DC_CD8_HIGH_LOW_CELLS.tsv.gz')
elig=[]
for r in regions.itertuples():
    for stratum in ['HIGH','LOW']:
        d=pair[(pair.region_id==r.region_id)&(pair.stratum==stratum)]
        ndc=int((d.source_native_cell_type==DC).sum());ncd8=int((d.source_native_cell_type==CD8).sum())
        elig.append(dict(patient_id=r.patient_id,region_id=r.region_id,state='Pre',stratum=stratum,n_DC=ndc,n_CD8=ncd8,n_union=len(d),count_evaluable=ndc>0 and ncd8>0,count_reason='EVALUABLE_COUNTS' if ndc>0 and ncd8>0 else 'MISSING_DC' if ndc==0 else 'MISSING_CD8',duplicate_XY_positions=int(d.duplicated(['x','y']).sum())))
elig=pd.DataFrame(elig);save(elig,'WU_AF35_DC_CD8_TOPOLOGY_ELIGIBILITY.tsv')
# Search bounded source-system inventories; ST-to-cell mappings and ST k=6 are not cell adjacency graphs.
graph_roots=[C,ROOT/'HCC_ICI_BIOLOGICAL_RESOLUTION_BR1_v1',ROOT/'HCC_ICI_project/03_processed_data/spatial_phase8',ROOT/'绘图视觉']
graph_candidates=[]
for folder in graph_roots:
    for p in folder.rglob('*'):
        if p.is_file() and any(k in p.name.lower() for k in ['adjacency','cell_graph','cell-graph','cell_edges','cell_neighbors']):graph_candidates.append(str(p))
assert not graph_candidates,'Inspect any frozen cell graph before protocol freeze'
write_json('QA/FROZEN_GRAPH_SEARCH.json',{'roots':[str(p) for p in graph_roots],'candidate_files':graph_candidates,'status':'NOT_RUN_NO_FROZEN_GRAPH','distinction':'Existing ST-to-cell mapping and directed k=6 ST graph are not CODEX cell-to-cell adjacency.'})
protocol=f'''# Wu AF35–DC–CD8 空间拓扑：结果前冻结方案 v1

冻结时间（UTC）：{now()}

## 授权与可行性

只读接续 Analysis C。其 54 个 manifest 条目已验证；所有重用文件先哈希再检查。151,302 个唯一映射细胞具有原生 CODEX X/Y、region_id、cell_id 和唯一冻结 AF35 stratum，逐细胞与 5,205 个区域点复核一致。未进行新配准、插值、最近邻映射、距离阈值选择或 ROI 选择。当前只检查计数与映射，尚未计算任何拓扑距离或 HIGH−LOW 拓扑结果。

患者：{', '.join(patients)}。全部 Pre。患者集保持不变；可估计性另行记录。

区域：
'''+'\n'.join(f'- {r.patient_id}: {r.region_id}' for r in regions.itertuples())+f'''

## 冻结定义

- 唯一细胞对：精确来源标签 `Dendritic cells` 与 `CD8 T cells`，不重标、不合并、不扫描其他细胞对。
- 坐标：Analysis C 单细胞表中的 CODEX 原始 X/Y，欧氏距离，单位为来源坐标单位；不将未核实的单位写成微米。各区域独立坐标，不跨区域连接。重复坐标保留；不移动、不抖动。
- HIGH/LOW：直接继承 Analysis C 的逐细胞 stratum。MIDDLE 不纳入本端点；保留既有母体清单。原 ST residual、type-7 四分位数和 mapping 均不修改。
- 最少计数：每个患者×区域×层至少 1 个 DC 和 1 个 CD8 才能数学定义 NND；不事后增加或降低阈值。稀疏层计数明示。Observed_NND 或 NullMedian_NND 非有限或 <=0 时该层不可估计，不补细胞、不加伪计数、不制造距离。
- 主 raw 指标：每个 DC 到同一区域同层最近 CD8 的欧氏距离，中位数，方向固定 DC→CD8。仅计算此方向，不增加反向端点。
- 空间零模型：固定 DC/CD8 占据位置的并集，只随机分配二元细胞类型标记；每次严格保持 DC 数、CD8 数、所有 XY 与 HIGH/LOW。每个可估计层 **2,000** 次置换，HIGH/LOW 相同。依据可行性计数最大并集 {int(elig.n_union.max())} 个位置，在查看距离结果前选定 2,000。
- seed 根值 **20260909**；NumPy PCG64。按患者数值、region_id 字典、LOW 再 HIGH 的固定层顺序，各层使用 SeedSequence([20260909, group_index])。每次从固定并集无放回抽取 n_DC 个位置标为 DC，其余 CD8；不移动坐标。完整置换索引及每次指标保存，逐次审核标签数量。
- 主分数 `log2(median(permutation_NND) / observed_NND)`；>0 表示比同计数、同并集占位的随机标记预期更近，<0 表示更远。0 为相等。未报告置换层面的 P 作为生物重复证据。
- 患者汇总：Analysis C 的细胞比例计数加权规则不适用于距离；按本请求等区域加权。分别取每位患者全部可估计 HIGH 区域分数的算术均值、全部可估计 LOW 区域分数的均值，随后 HIGH−LOW。两个层使用的区域集合可以不同，但明确列出；不跨区域计算邻居，不按细胞计数加权。raw 距离用相同等区域规则。不可估计层保留原因。
- >=8 位同时具备 HIGH 与 LOW 的患者才执行主推断；不足即停止推断，状态 NOT_EVALUABLE_LOW_PATIENT_COVERAGE，不放宽条件。
- 主推断：所有 2^N 患者符号配置，双侧检验均值绝对值，使用 >= observed−1e-12 数值容差，包含全部零效应患者。该符号交换零假设为统计假设，非随机化因果设计。
- 辅助区间：患者等权 10,000 次 bootstrap，PCG64 seed **20260910**，mean 的 percentile 2.5%/97.5% 区间，NumPy linear/type-7 分位数。CI 辅助，exact P 主判定。保存 bootstrap 患者索引以独立复现。
- 成功规则：N>=8、正向数 >=ceil(0.8*N)、均值>0、exact P<0.05；此 DELTA 已为丰度保留零模型调整效应，满足请求 E。N=10 时需至少8正向。正均值但不满足完整成功标准记 DIRECTIONALLY_SUPPORTIVE_NOT_SIGNIFICANT（若仅方向一致性未过，另写确切原因）；均值<=0记 NOT_SUPPORTED。不用 raw 结果覆盖主结论。
- 次要分析：患者 raw DC→CD8 NND，方向 HIGH<LOW，描述性。不执行新 raw 假设检验。未找到既有冻结 CODEX cell-to-cell adjacency graph；DC_CD8_ADJACENCY_SECONDARY=NOT_RUN_NO_FROZEN_GRAPH。现有 ST 图和 ST→cell mapping 不能替代细胞邻接图。
- 结果不是 Analysis C 的独立队列复制；仅同患者系统的预指定拓扑延伸。可支持相邻切片空间关联，不能证明物理相互作用、抗原呈递、招募、通信、因果机制或新生态位。

## 验证与停止

独立实现：SciPy cKDTree 主 NND，R 直接二维欧氏距离复算 observed；R 等区域患者 DELTA、exact sign-flip 与保存抽样索引的 bootstrap 独立复核。每次置换检查标记计数与坐标不变，保存索引可重算。图的患者点与表逐值绑定。保护旧 AF35、Analysis C、A/B 等补充分析、主图、图注与手稿，前后 SHA-256 比对，不写入旧目录。

完成后 STOP，按实际证据返回预定 gate，不启动 Analysis E，不更换指标、细胞对、距离阈值或置换次数，不修改 Figure 5、Results、图注或手稿。
'''
(OUT/'WU_AF35_DC_CD8_TOPOLOGY_PROTOCOL_FREEZE_v1.md').write_text(protocol,encoding='utf-8')
bound={r['path']:r['sha256'] for r in binding}
bound[str(OUT/'data/FROZEN_DC_CD8_HIGH_LOW_CELLS.tsv.gz')]=sha(OUT/'data/FROZEN_DC_CD8_HIGH_LOW_CELLS.tsv.gz')
lock={'freeze_time_utc':now(),'protocol_sha256':sha(OUT/'WU_AF35_DC_CD8_TOPOLOGY_PROTOCOL_FREEZE_v1.md'),'bound_inputs':bound,'patients':patients,'regions':regions[['patient_id','region_id']].to_dict('records'),'permutations':2000,'permutation_seed':20260909,'bootstrap_n':10000,'bootstrap_seed':20260910,'primary_direction':'DC_to_CD8','outcome_analysis_started':False,'feasibility_gate':'PASS'}
write_json('WU_AF35_DC_CD8_TOPOLOGY_FREEZE_LOCK.json',lock)
write_json('QA/FEASIBILITY_RESULT.json',{'gate':'PASS','n_C_cells':len(edges),'n_pair_high_low_cells':len(pair),'n_region_strata':len(elig),'n_count_evaluable':int(elig.count_evaluable.sum()),'max_union':int(elig.n_union.max()),'protected_files':len(records),'topology_calculated_before_freeze':False})
print(elig.to_string(index=False),flush=True)
print('FEASIBILITY_PASS; PROTOCOL_FROZEN; NO_TOPOLOGY_CALCULATED',now(),flush=True)
