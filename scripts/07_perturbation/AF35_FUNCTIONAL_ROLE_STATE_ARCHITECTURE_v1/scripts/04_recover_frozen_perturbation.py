import importlib,sys,io
from collections import OrderedDict
import numpy as np
f=importlib.import_module('00_initialize');OUT=f.OUT
sys.path.insert(0,str(f.P/'10_software/phase3_py'));import h5py
genes=f.read(OUT/'inputs/AF35_ORDERED_GENES.tsv').gene.tolist()
base=f.P/'10_intermediate_files/phase17a_public_perturbation'
saved=f.jl(base/'PHASE17A_ANALYSIS_RESULTS.json');src=saved['gse314342']['remote_sources'][0]
cache=base/'GSE314342/combined_h5ad_range_cache'
class CacheOnly(io.RawIOBase):
 def __init__(self):self.position=0;self.size=src['size_bytes'];self.blocks=OrderedDict();self.used=set()
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.position
 def seek(self,o,w=0):self.position=o if w==0 else (self.position+o if w==1 else self.size+o);return self.position
 def read(self,size=-1):
  end=self.size if size<0 else min(self.size,self.position+size);chunks=[]
  while self.position<end:
   n=self.position//(8*1024*1024);offset=self.position%(8*1024*1024)
   if n not in self.blocks:
    p=cache/f'{src["etag"]}_{n:08d}.bin'
    if not p.is_file():raise FileNotFoundError('FROZEN_CACHE_BLOCK_MISSING:'+str(p))
    self.blocks[n]=p.read_bytes();self.used.add(p)
    if len(self.blocks)>8:self.blocks.popitem(last=False)
   block=self.blocks[n];take=min(end-self.position,len(block)-offset);assert take>0
   chunks.append(block[offset:offset+take]);self.position+=take
  return b''.join(chunks)
 def readinto(self,b):v=self.read(len(b));b[:len(v)]=v;return len(v)
def decode(v):return np.array([s.decode() if isinstance(s,bytes) else str(s) for s in v])
def field(h,n):
 x=h[n]
 return decode(x[:]) if isinstance(x,h5py.Dataset) else decode(x['categories'][:])[x['codes'][:]]
status='NOT_AVAILABLE_AS_FROZEN_GENE_LEVEL_ENDPOINT';detail='';rows=[];checks=[];reader=CacheOnly()
try:
 with h5py.File(reader,'r') as h:
  tg=field(h,'obs/target_contrast_gene_name');cond=field(h,'obs/culture_condition');gn=decode(h['var/gene_name'][:]);target=np.flatnonzero(tg=='PTPN22')
  assert len(target)==3 and all(np.sum(gn==g)==1 for g in genes)
  mapping={'Rest':'PTPN22_CRISPRI_VS_NTC_REST','Stim8hr':'PTPN22_CRISPRI_VS_NTC_STIM8HR','Stim48hr':'PTPN22_CRISPRI_VS_NTC_STIM48HR'}
  for i in target:
   for j,g in enumerate(genes):
    col=np.flatnonzero(gn==g).item()
    lfc=float(h['layers/log_fc'][i,col]);z=float(h['layers/zscore'][i,col]);padj=float(h['layers/adj_p_value'][i,col])
    rows.append(dict(gene_rank=j+1,gene=g,dataset='GSE314342',contrast=mapping[cond[i]],condition=cond[i],frozen_author_log_fc=lfc,direction='POSITIVE' if lfc>0 else ('NEGATIVE' if lfc<0 else 'ZERO'),author_zscore=z,author_adjusted_p=padj,donor_n=4,inference_unit='AUTHOR_DONOR_AWARE_POOLED_DE_MODEL',source_url=src['url'],source_etag=src['etag'],source_type='EXISTING_LOCAL_FROZEN_AUTHOR_DE_CACHE',new_endpoint_fitted=False))
 for condition,contrast in mapping.items():
  r=[r for r in rows if r['condition']==condition];old=next(r for r in saved['gse314342']['program_effects'] if r['contrast_id']==contrast and r['program']=='PTPN22_activation_feedback_35')
  values=np.array([r['frozen_author_log_fc'] for r in r]);zs=np.array([r['author_zscore'] for r in r]);ps=np.array([r['author_adjusted_p'] for r in r])
  for key,value in [('mean_author_log_fc',values.mean()),('median_author_log_fc',np.median(values)),('mean_author_zscore',zs.mean()),('same_sign_gene_fraction',np.mean(np.sign(values)==np.sign(values.mean()))),('adjusted_p_lt_0_05_gene_n',int((ps<.05).sum())),('present_genes',len(values))]:
   err=abs(value-old[key]);checks.append(dict(contrast=contrast,metric=key,recovered=value,frozen_reference=old[key],absolute_error=err,status='PASS' if err<1e-10 else 'FAIL'))
 assert all(r['status']=='PASS' for r in checks)
 status='PASS';f.save(rows,'AF35_PTPN22_PERTURBATION_GENE_EFFECTS.tsv')
 detail='从既有本地作者DE缓存恢复35基因×3个原对比；不请求网络字节，不拟合DE，不改contrast。全部原AF35程序均值、中位数、z、同向比例、调整P计数和覆盖复现。'
except (FileNotFoundError,OSError,AssertionError) as e:
 detail='来源恢复不足，未补救：'+repr(e)
f.save(checks if checks else [dict(status='NOT_RUN',detail=detail)],'QA/PERTURBATION_FROZEN_AGGREGATE_VERIFICATION.tsv')
files=[base/'PHASE17A_ANALYSIS_RESULTS.json',f.P/'05_results/tables/phase17a_public_perturbation/PHASE17A_GSE314342_PROGRAM_EFFECTS.tsv',f.P/'04_scripts/phase17_dual_strengthening/02_analyze_phase17a.py',f.P/'05_results/phase17a_public_perturbation/PHASE17A_ANALYSIS_METHOD_FREEZE.md']+sorted(reader.used)
before=f.read(OUT/'QA/PROTECTED_BEFORE.tsv').set_index('path')
bound=[]
for p in files:
 actual=f.sha(p);assert actual==before.loc[str(p),'sha256']
 bound.append(dict(path=str(p),sha256=actual,source_binding='PRE_I_PROTECTED_BASELINE_AND_ORIGINAL_ETAG'))
f.save(bound,'QA/PERTURBATION_SOURCE_BINDINGS.tsv')
f.js('QA/PERTURBATION_RECOVERY_LOCK.json',dict(utc=f.now(),gate=status,details=detail,source_etag=src['etag'],cached_blocks_used=len(reader.used),remote_bytes_requested=0,new_DE_models=0))
f.write('AF35_PTPN22_PERTURBATION_RECOVERY_AUDIT.md',f'''# 既有Figure4扰动基因端点恢复

AF35_PTPN22_PERTURBATION_RECOVERY_GATE = {status}

{detail}

范围为当前F4 ledger的GSE265961机制及Phase17扰动系统。GSE265961现有结果是condition/readout及程序汇总；GSE264373也是程序得分及重复层效应，本次不从原counts重建35基因DE。
GSE314342原Phase17A方法明确用作者donor-aware pooled DE的layers/log_fc、zscore和adj_p_value汇总预冻结AF35。原程序脚本在既有3个PTPN22行读取这些基因统计，因此此处是缓存中已有作者DE端点的提取，不是新拟合端点。缓存只读，size/etag来自已冻结PHASE17A_ANALYSIS_RESULTS.json；本轮逐块哈希与初始保护基线比较，原AF35聚合值再绑定数值身份。缓存文件原Phase17清单未逐块给SHA，故历史绑定同时依靠记录的etag和原冻结聚合复现；不夸大为全远程对象逐字节历史校验。
只用Rest/Stim8hr/Stim48hr原对比，不补缺失供者子集，不将时间条件作为AF35层级演进轨迹。作者adjusted_p与z为原不确定性/检验度量，不是本次新检验或效应CI。保留所有35基因的正负及不显著结果，donor_n=4不是细胞数。最多称perturbation-linked support，不能推出PTPN22直接控制AF35级联。
''')
print('I perturbation recovery:',status,detail,flush=True)
