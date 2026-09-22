from pathlib import Path
import hashlib, json, math, datetime, shutil
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent
P = ROOT / 'HCC_ICI_project'
def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1048576), b''): h.update(b)
    return h.hexdigest()
def tsv(rows, name): pd.DataFrame(rows).to_csv(OUT / name, sep='\t', index=False)
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()

if __name__ == '__main__':
    assert not (OUT / 'FREEZE_LOCK.json').exists(), 'Already frozen; do not overwrite'
    files = [p for p in ROOT.rglob('*') if p.is_file() and OUT not in p.parents]
    tsv([dict(path=str(p.relative_to(ROOT)), bytes=p.stat().st_size, mtime_ns=p.stat().st_mtime_ns) for p in files], 'PROTECTED_EXISTING_FILE_INVENTORY.tsv')
    sources = [
      P/'04_scripts/phase17_dual_strengthening/03_analyze_phase17b_outcome_blind.py',
      P/'04_scripts/phase15c_public_benchmark/01_run_phase15c_benchmark.R',
      P/'01_metadata/spatial_phase8/PTPN22_SPATIAL_PROGRAM_FROZEN_PROVENANCE.tsv',
      P/'01_metadata/clinical_translation/PTPN22_CLINICAL_AXES_FROZEN.tsv',
      P/'03_processed_data/pretreatment_clinical_translation/ERP117672_gene_expression_response_blind.rds',
      P/'03_processed_data/phase15b_public_replication/GSE302495_RESPONSE_BLIND_MATRIX_AND_SCORES.rds',
      P/'01_metadata/pretreatment_clinical_translation/ERP117672_PATIENT_RESPONSE_MAP.tsv',
      P/'01_metadata/phase15b_public_replication/GSE302495_PATIENT_RESPONSE_MAP.tsv',
      P/'01_metadata/phase17_dual_strengthening/PHASE17_FILE_MANIFEST_SHA256.tsv',
      ROOT/'HCC_ICI_LATEST_PROJECT_CHECKPOINT.md',
      ROOT/'HCC_ICI_checkpoint_2026-08-31_v19_FINAL_PREWRITING_STRENGTHENING_FREEZE.md',
      ROOT/'HCC_ICI_checkpoint_2026-08-31_v20_CONTROLLED_MANUSCRIPT_CONSTRUCTION_FREEZE.md']
    sources += list((P/'05_results/tables/phase17b_af35_portability').glob('*.tsv'))
    sources += list((P/'05_results/phase17b_af35_portability').glob('*.md'))
    sources += [P/f'05_results/tables/phase15c_public_benchmark/{c}_UNIVARIATE_BENCHMARK_EFFECTS.tsv' for c in ['ERP117672','GSE302495']]
    sources += [p for p in files if any(x in p.name.upper() for x in ['MANUSCRIPT','MAIN_FIGURE','LEGEND','RESULTS_DRAFT']) and p.suffix.lower() in ['.md','.tsv','.docx','.svg','.pdf','.png']]
    sources = sorted(set(sources))
    tsv([dict(path=str(p.relative_to(ROOT)), bytes=p.stat().st_size, sha256=sha(p)) for p in sources], 'PROTECTED_INPUT_SHA256_BEFORE.tsv')
    expected = {'ERP117672_gene_expression_response_blind.rds':'ec4c1fd4c25d04e81167766c252da72047256a0352cd416920aba74847cc63f6', 'GSE302495_RESPONSE_BLIND_MATRIX_AND_SCORES.rds':'b138295da8af3f7f8965dc68ff429aa10cf63e90f3e0ce286a1185e299322d1e','ERP117672_PATIENT_RESPONSE_MAP.tsv':'4f28ae48fabf509f1248ac808d2c07b1a6beaf37b4f8f14fbeb999431f7082b2','GSE302495_PATIENT_RESPONSE_MAP.tsv':'fd2d8101719bd178e0b28680645e86319a30e485f210dd8f27fae04dc1c83b8d'}
    bindings = []
    for p in sources:
        if p.name in expected:
            actual = sha(p)
            assert actual == expected[p.name], f'SOURCE_AUTHORITY_CONFLICT: {p}'
            bindings.append(dict(source=str(p.relative_to(ROOT)), expected_sha256=expected[p.name], actual_sha256=actual, status='PASS', outcome_content_read='NO'))
    tsv(bindings, 'AF35_SOURCE_BINDING_AUDIT.tsv')
    g = pd.read_csv(P/'01_metadata/spatial_phase8/PTPN22_SPATIAL_PROGRAM_FROZEN_PROVENANCE.tsv', sep='\t').sort_values('gene_rank').gene.tolist()
    axes = pd.read_csv(P/'01_metadata/clinical_translation/PTPN22_CLINICAL_AXES_FROZEN.tsv', sep='\t')
    assert g == axes.loc[axes.component.eq('PTPN22_activation_feedback_35'),'exact_features'].iloc[0].split(';')
    assert len(g) == len(set(g)) == 35 and 'PTPN22' not in g
    tsv([dict(gene_rank=i+1,gene=x) for i,x in enumerate(g)], 'AF35_ORDERED_GENES.tsv')
    tsv([dict(subset_id=f'SINGLE_{i+1:02}', omitted_gene=x, removed_gene_n=1, retained_gene_n=34, retained_genes=';'.join(y for y in g if y!=x)) for i,x in enumerate(g)], 'AF35_SINGLE_GENE_SUBSET_MANIFEST.tsv')
    rows=[]
    for ci,c in enumerate(['ERP117672','GSE302495'],1):
        for pct,k in [(10,4),(20,7),(30,11),(40,14),(50,18)]:
            seed=20260830+ci*100000+pct*10
            rng=np.random.default_rng(seed)
            for draw in range(1,1001):
                removed=set(rng.choice(g,size=k,replace=False).tolist())
                rows.append(dict(cohort=c,subset_id=f'{c}_D{pct}_{draw:04}',draw=draw,dropout_percent=pct,removed_gene_n=k,retained_gene_n=35-k,retained_fraction=(35-k)/35,actual_removed_fraction=k/35,rng_seed=seed,omitted_genes=';'.join(x for x in g if x in removed),retained_genes=';'.join(x for x in g if x not in removed),origin='RECOVERED_PHASE17B_SEED' if pct<=30 else 'RESPONSE_BLIND_STRESS_EXTENSION'))
    tsv(rows,'AF35_RANDOM_DROPOUT_SUBSET_MANIFEST.tsv')
    (OUT/'AF35_COMPONENT_ABLATION_PROTOCOL_FREEZE_v1.md').write_text('''# AF35 组分消融方案冻结 v1

本轮用户附件为新增分析授权；仅在此独立目录执行，旧稿件、主图、Results 和既有分析不作修改。PTPN22 为 BIOLOGICAL_ANCHOR；AF35 为固定无符号、等权重、排除 PTPN22 的 35 基因 CLINICAL_STATE_UNIT。子集仅为组分扰动，不能命名为新签名、挑选最优子集或赋予因果重要性。

## 结局防火墙与来源
在接入患者结局之前生成、哈希全部子集和本方案；已知用户提供的总体效应和事件数，不声称分析者对既往结果不知情。恢复阶段仅读取无结局表达对象、基因定义、既有技术统计，结局映射文件只作二进制哈希校验，不解析内容。每一阶段保存 UTC 时间和哈希。

## 评分与分析人群
ERP117672：冻结 log2(TPM+1) 对象全部 40 例上逐基因均值中心化、样本标准差 (ddof=1) 标准化，然后取保留基因等权平均并在全部 40 例重标准化。关联计算才限制为原映射中 primary_included 的 35 例 (6/29)。GSE302495：冻结 log2(CPM+1) 对象全部 38 例上采用同样评分 (12/26)。删除基因不重估其他基因的标准化参数；无缺失插补。完整 AF35 与原始冻结效应必须一致 (绝对误差 <1e-10)。

## 子集抽样
35 个单基因删除，保留原始顺序。10/20/30% 精确恢复 Phase17B 的 NumPy default_rng/PCG64、队列序号 ERP=1/GSE=2、seed=20260830+队列序号*100000+比例*1000，删除数 4/7/11。原实现硬编码删除数，没有通用舍入函数；其数值符合 ceiling(35*fraction)，本次扩展明确冻结此规则，40/50% 删除 14/18。每档每队列 1000 次，子集内部无放回，不同抽样之间可重复且不去重。不同档独立抽样，非嵌套删除路径；两队列沿用不同子集，不作配对子集比较。10/20/30% 必须用既有 rho 摘要和所有 35 个 LOO rho 验证重现 (<1e-12)，验证仅为恢复 QA，不作为新技术发现。

## 效应与统计范围
Hedges g = (1-3/(4*(n-2)-1))*(mean_R-mean_NR)/pooled_sample_SD，与冻结实现一致。delta_g 相对同队列完整 AF35。rho 主表使用该队列关联分析患者；恢复 QA 的 ERP rho 使用全部 40 例，与既有技术统计保持对应。方向严格按 >0/<0，零值记 ZERO，解释为非正。分位数采用 NumPy linear。随机子集 2.5/25/50/75/97.5 百分位是 COMPONENT-PERTURBATION 分布；不是患者抽样置信区间或生物重复。无多重检验、p 值筛选、ROC/AUC、阈值优化、分类器、权重学习、跨队列合并或外部队列扩展。

## 预设解释规则
以下为本轮描述性审查规则，不是已验证生物学或临床阈值，按较不利队列形成模块标签。SINGLE_GENE_DEPENDENCE：任一删除 g<=0 或 |delta_g|>=0.30 为 SINGLE_GENE_SENSITIVE；否则所有删除 |delta_g|<=0.15 且 rho>=0.95 为 DISTRIBUTED；其余 MODERATE_COMPONENT_SENSITIVITY。
RANDOM_DROPOUT_BEHAVIOR：任一队列在任一主档 (10/20/30%) 正向比例<0.80 或中位 g<完整 g 的 50% 为 EARLY_COLLAPSE；否则所有队列全部五档正向比例>=0.95、中位 g>=完整 g 的 75% 且中位 rho>=0.90 为 GRACEFUL_DEGRADATION；其余 INTERMEDIATE。类别名称不预设下降必须单调，需另报实际趋势。
CROSS_COHORT_DIRECTION_STABILITY：两个队列均 35/35 单删正向且全部主档正向比例>=0.95 为 STRONG；若任一队列所有单删及主档比例均>=0.80 为 PARTIAL；否则 WEAK。另报同一基因删除在任一/两个队列逆向的数量。

## 可选锚点/读出条件敏感性（本轮计划运行）
复用 R logistf 冻结框架：response ~ standardized PTPN22 + standardized AF35。在各自最终可评估患者中两预测量再按样本 SD 标准化，Firth=TRUE，pl=TRUE；报告两系数、每 SD OR、95% profile penalized-likelihood 区间、事件数及事件/参数。无模型优越性检验或预测指标。任一拟合失败、非有限系数/区间、max(abs(conv))>=1e-3 或 profile 迭代达到上限则 NOT_INTERPRETABLE；成功但事件/参数<10 时 PASS_WITH_LIMITS。两队列 AF35 系数均>0 为 POSITIVE_BOTH，否则 MIXED；未通过数值条件为 NOT_INTERPRETABLE。正向不等于区间排除零效应；该结果不支持机制独立或预测改进。

## 终止
全部有效性、来源、防火墙、评分、哈希及图件核验完成后，最终门为 PASS_FOR_AUTHOR_SCIENTIFIC_REVIEW，随即停止。本方案和子集一经冻结不得根据新结果调整。
''',encoding='utf-8')
    pass  # Non-computational private authorization copy excluded.
    frozen=['AF35_ORDERED_GENES.tsv','AF35_SINGLE_GENE_SUBSET_MANIFEST.tsv','AF35_RANDOM_DROPOUT_SUBSET_MANIFEST.tsv','AF35_COMPONENT_ABLATION_PROTOCOL_FREEZE_v1.md']
    hashes={n:sha(OUT/n) for n in frozen}
    tsv([dict(path=n,sha256=h) for n,h in hashes.items()],'AF35_SUBSET_MANIFEST_SHA256.tsv')
    (OUT/'FREEZE_LOCK.json').write_text(json.dumps(dict(time_utc=now(),outcome_labels_loaded=False,hashes=hashes,numpy_version=np.__version__),indent=2),encoding='utf-8')
    print(json.dumps(dict(stage='MANIFEST_FROZEN',source_hashes=len(sources),existing_files=len(files),subsets=len(rows)+35)))
