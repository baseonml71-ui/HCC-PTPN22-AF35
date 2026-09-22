from pathlib import Path
import csv, gzip, hashlib, importlib.util, json, shutil, sys, zipfile
import numpy as np
import pandas as pd
sys.dont_write_bytecode = True

OUT = Path(__file__).resolve().parents[1]
ROOT = OUT.parent
PROJECT = ROOT / 'HCC_ICI_project'
def sha(p):
    with open(p, 'rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()
def save(rows, name):
    pd.DataFrame(rows).to_csv(OUT / name, sep='\t', index=False, na_rep='NA')

def main():
    if (OUT / 'PROTECTED_INPUT_BEFORE.tsv').exists():
        raise RuntimeError('Authority baseline already exists; do not overwrite the freeze')
    for d in ['authority', 'reproduction', 'qa', 'figures', 'scripts']:
        (OUT / d).mkdir(exist_ok=True)
    packet = ROOT / '补充分析‘/HCC_ICI_D_CLONE_STATE_DECOMPOSITION_PACKET_2026-09-08.zip'
    with zipfile.ZipFile(packet) as z:
        for e in z.infolist():
            target = (OUT / 'authority' / e.filename).resolve()
            assert target.is_relative_to(OUT / 'authority') and not e.is_dir()
            target.write_bytes(z.read(e))
    source_paths = [
        '04_scripts/phase19/02_phase19a_patient_analysis.py',
        '04_scripts/phase19/01_phase19a_response_blind_projection.py',
        '04_scripts/phase3/03_independent_scrna_sctcr_clone_state.R',
        '04_scripts/phase_temporal/01_gse235863_baseline_future_clone_fate.R',
        '03_processed_data/phase3/GSE235863_BARCODE_CLONOTYPE_STATE_PTPN22.tsv.gz',
        '03_processed_data/phase19a_af35_clonotype/GSE235863_AF35_RESPONSE_BLIND_PROJECTION.npz',
        '01_metadata/phase19a_af35_clonotype/PHASE19A_AF35_RESPONSE_BLIND_PROJECTION_FREEZE.json',
        '01_metadata/phase19a_af35_clonotype/PHASE19A_AF35_SINGLE_CELL_PROJECTION_CONTRACT.md',
        '01_metadata/phase19_final_prewriting/PHASE19_FILE_MANIFEST_SHA256.tsv',
        '01_metadata/phase_temporal/GSE235863_V6_TEMPORAL_FREEZE_MANIFEST_SHA256.tsv',
        '01_metadata/final_harmonization/GSE235863_TEMPORAL_COMPARTMENT_AUDIT.tsv',
        '05_results/temporal/GSE235863_LONGITUDINAL_CLONOTYPE_ID_AUDIT.md',
        '05_results/tables/temporal/GSE235863_BASELINE_CLONE_FATE_TABLE.tsv',
        '05_results/tables/phase3/GSE235863_WITHIN_PATIENT_CLONAL_ASSOCIATIONS.tsv',
        '05_results/tables/phase3/GSE235863_CLONE_PATIENT_SUMMARIES.tsv',
        '05_results/tables/phase19a_af35_clonotype/AF35_CLONOTYPE_TEMPORAL_SENSITIVITY.tsv',
        '05_results/tables/phase19a_af35_clonotype/AF35_EXPANDED_SINGLETON_PATIENT_EFFECTS.tsv',
        '05_results/tables/phase19a_af35_clonotype/AF35_PERSISTENT_NONPERSISTENT_PATIENT_EFFECTS.tsv',
        '05_results/tables/phase19a_af35_clonotype/AF35_BASELINE_TO_FUTURE_OBSERVED_PATIENT_EFFECTS.tsv',
        '10_intermediate_files/phase19/phase19a_table_data.json',
        '10_intermediate_files/phase19/phase19a_summary.json',
    ]
    protected = {PROJECT / p for p in source_paths} | {packet}
    # Protect existing manuscript text and Figure 3 files by content, without claiming a whole-workspace lock.
    protected |= {p for p in ROOT.rglob('*.md') if not p.is_relative_to(OUT) and
                  ('MANUSCRIPT' in p.name.upper() or 'CHECKPOINT' in p.name.upper() or 'FIGURE3' in str(p).upper() or 'F3_' in p.name)}
    for d in (ROOT / '绘图视觉').iterdir():
        if d.is_dir() and any(t in d.name.upper() for t in ['FIGURE3', '_F3', 'V2_1H']):
            protected |= {p for p in d.rglob('*') if p.is_file() and p.suffix.lower() in ['.svg','.pdf','.png','.md','.tsv','.r','.py','.json']}
    records = []
    for p in sorted(protected):
        assert p.exists(), p
        records.append(dict(path=str(p.relative_to(ROOT)), bytes=p.stat().st_size, sha256=sha(p)))
    save(records, 'PROTECTED_INPUT_BEFORE.tsv')
    print('Protected hashes:', len(records), flush=True)
    for p in source_paths:
        src = PROJECT / p
        dst = OUT / 'authority' / p
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    freeze = json.loads((PROJECT / source_paths[6]).read_text(encoding='utf-8'))
    assert sha(PROJECT / source_paths[4]) == freeze['mapping_sha256']
    assert sha(PROJECT / source_paths[5]) == freeze['npz_sha256']
    assert freeze['af35_gene_count'] == freeze['usable_gene_count'] == 35
    assert freeze['outcome_joined'] is False and freeze['ptpn22_excluded'] is True
    manifest = pd.read_csv(PROJECT / source_paths[8], sep='\t')
    binding = []
    for rel in source_paths:
        row = manifest[manifest.relative_path == 'HCC_ICI_project/' + rel]
        if len(row):
            actual = sha(PROJECT / rel)
            expected = row.iloc[0].sha256
            binding.append(dict(path=rel, expected_sha256=expected, actual_sha256=actual, status='PASS' if expected == actual else 'FAIL'))
    save(binding, 'HISTORICAL_SOURCE_HASH_VERIFICATION.tsv')
    assert all(r['status']=='PASS' for r in binding)
    (OUT / 'CLONE_STATE_DECOMPOSITION_PROTOCOL.md').write_text('''# Clone-State Decomposition 执行方案

授权来源为 2026-09-08 Analysis D 包。9 月 6 日交接文件只作背景，不执行旧 Figure 3 任务。
仅使用冻结本地输入，不读取疗效标签进行分组、筛选或调参。原脚本复现保留其原有输出但不开展新的疗效分析。
患者为生物学重复。扩增使用全部 58,872 个细胞、作者跨样本 clone.size >= 2；持续检出使用 persistence_eligible 的精确患者-组织集合；后来检出使用既定基线外周血集合。
校正算子为每位患者内 sub_cluster 分层差、min(n_positive,n_reference) 加权；不做新回归、残差化或重标准化。没有两组同时存在的层不进入校正，输出纳入/排除细胞数与层明细。
PTPN22 时间比较本轮采用与 AF35 相同的基线细胞加权算子，以保证 raw/adjusted 可比；这与 v6 的每克隆等权、患者中位数估计不同，单列并明确标记，绝不替换 v6 冻结结果。
原始脚本在独立 reproduction 目录完整重放，保留种子 20260831 和 RNG 调用顺序以核验冻结区间。新增比较使用相同 10,000 次患者 bootstrap 百分位法，种子 20260908；差异仅来自预先固定的新随机流，历史区间另列。
衰减 difference=raw-adjusted，ratio=adjusted/raw，fraction=1-ratio，纯描述、不解释为中介。输出同一可比较状态支持集的 raw 对照，提示校正同时改变权重/支持集。
重叠先定义资格：扩增×持续检出按患者-组织克隆行；扩增×后来检出按基线外周血克隆行；持续检出×后来检出及三重重叠仅在相同基线外周血集合中计算。不同组织不合并成统一克隆类别。患者无条件汇总不作为新生物学推断。
配对 CDR3 敏感性在患者内重建全细胞克隆大小及同组织 pre/post 检出，保留主分析细胞资格；作者 clone.id 始终为主。
肿瘤/血液可行性先按同患者、同 pre/post、GEO 精确样本标题验证。采用匹配时点两个组织的克隆并集；同一患者多个时点不增加 n。技术 PASS 要求至少 3 位患者有来源有效的匹配时点，且两组织内均存在共检出与局限检出可比较细胞；不依赖表达方向。先输出可行性与 clone-size 审计，再运行 Supplement-only 患者内、组织内比较。没有独立深度控制时不作迁移或交通推断。
主科学图保持独立尺度，无轨迹箭头。全部产物仅写本目录，完成后 STOP_FOR_AUTHOR_SCIENTIFIC_REVIEW。
''', encoding='utf-8')
    spec = importlib.util.spec_from_file_location('frozen_phase19', PROJECT / source_paths[0])
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.WORK_DIR = OUT / 'reproduction'
    mod.TABLE_JSON = mod.WORK_DIR / 'phase19a_table_data.json'
    mod.SUMMARY_JSON = mod.WORK_DIR / 'phase19a_summary.json'
    mod.LOG_FILE = mod.WORK_DIR / '02_phase19a_patient_analysis.log'
    mod.main()
    original = json.loads((PROJECT / source_paths[-2]).read_text(encoding='utf-8'))
    replay = json.loads(mod.TABLE_JSON.read_text(encoding='utf-8'))
    checks = []
    def compare(a,b,path='root'):
        if isinstance(a,dict):
            assert a.keys()==b.keys(), path
            for k in a: compare(a[k],b[k],path+'.'+k)
        elif isinstance(a,list):
            assert len(a)==len(b),path
            for i,(x,y) in enumerate(zip(a,b)): compare(x,y,path+f'[{i}]')
        elif isinstance(a,(float,int)) and not isinstance(a,bool):
            ok = (np.isnan(a) and np.isnan(b)) or np.isclose(a,b,rtol=0,atol=1e-12)
            checks.append(dict(field=path, original=a, replay=b, status='PASS' if ok else 'FAIL'))
        else: assert a==b, (path,a,b)
    compare(original,replay)
    save(checks, 'FROZEN_PHASE19_NUMERIC_REPRODUCTION.tsv')
    assert all(r['status']=='PASS' for r in checks)
    print('Authority replay verified:', len(checks), flush=True)

if __name__ == '__main__':
    main()
