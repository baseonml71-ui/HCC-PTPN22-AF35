from common import *
from pandas.io.formats import csvs
assert not (Q/'PROJECTION_COMPLETION.json').exists()
gate=jl(Q/'SOURCE_GATE_LOCK.json');assert gate['background']=='PASS'
allowed={str((OUT/(P+'_PERTURBATION_BACKGROUND_UNIVERSE.tsv')).resolve()),str((A/'AF35_ORDERED_GENES.tsv').resolve()),str((OUT/(P+'_PERTURBATION_ORIENTATION_FREEZE.tsv')).resolve())}
access=[]
def guard(event,args):
 if event=='open' and isinstance(args[0],(str,bytes,os.PathLike)):
  path=Path(os.fsdecode(args[0])).resolve();mode=str(args[1]);flags=args[2] or 0;writing=any(c in mode for c in 'wax+') or bool(flags&(os.O_WRONLY|os.O_RDWR|os.O_CREAT))
  ok=path.is_relative_to(OUT) if writing else str(path) in allowed
  access.append(dict(path=str(path),writing=writing,allowed=ok))
  if not ok:raise PermissionError('N_CLINICAL_FIREWALL_DENIED: '+str(path))
 if event.startswith(('socket.','subprocess.')) or event in ['os.system','os.spawn']:raise PermissionError('N_NETWORK_PROCESS_DENIED')
sys.addaudithook(guard)
for p in [I/'README.md',ROOT/'AF35_LONGITUDINAL_CLONAL_STATE_DYNAMICS_M2_CLEANROOM_v1/03_phaseB_outcome_inputs/RESPONSE.tsv']:
 try:open(p,'rb');raise AssertionError('Guard failed')
 except PermissionError:pass
bg=read(OUT/(P+'_PERTURBATION_BACKGROUND_UNIVERSE.tsv'));og=read(A/'AF35_ORDERED_GENES.tsv');genes=og.gene.tolist();ori=read(OUT/(P+'_PERTURBATION_ORIENTATION_FREEZE.tsv'));systems=ori.system.tolist()
assert len(genes)==35 and len(set(genes))==35 and 'PTPN22' not in genes
bg=bg[bg.eligible].copy();bg['oriented_effect']=bg.native_effect*bg.orientation_factor
def rank(v):return 2*(rankdata(v,method='average')-.5)/len(v)-1
def reference_summary(v,obs,identity):
 below=int(np.sum(v<obs-1e-12));equal=int(np.sum(abs(v-obs)<=1e-12))
 return dict(reference=identity,observed=obs,random_median=float(np.median(v)),random_q025=float(np.quantile(v,.025)),random_q975=float(np.quantile(v,.975)),random_n=10000,count_ge=int(np.sum(v>=obs-1e-12)),count_le=int(np.sum(v<=obs+1e-12)),count_equal=equal,reference_percentile=100*(below+.5*equal)/10000,formal_p_value_reported=False)
allranks=[];program=[];refall=[];refsummary=[];pools=[]
for i,system in enumerate(systems):
 d=bg[bg.system.eq(system)].sort_values('gene').copy();d['rank_score']=rank(d.oriented_effect.to_numpy());d['average_rank']=rankdata(d.oriented_effect.to_numpy(),method='average');d['percentile_rank']=(d.average_rank-.5)/len(d)
 af=d.set_index('gene').loc[genes].copy();af['gene_rank']=np.arange(1,36);allranks.append(af.reset_index());obs=float(af.rank_score.mean())
 program.append(dict(system=system,condition=ori.condition.iloc[i],eligible_genes=len(d),AF35_genes=35,AF35_PROGRAM_SHIFT=obs,mean_raw_oriented_effect=float(af.oriented_effect.mean()),median_raw_oriented_effect=float(af.oriented_effect.median()),positive_gene_n=int((af.oriented_effect>0).sum()),negative_gene_n=int((af.oriented_effect<0).sum()),zero_gene_n=int((af.oriented_effect==0).sum()),positive_fraction=float((af.oriented_effect>0).mean()),negative_fraction=float((af.oriented_effect<0).mean())))
 pool=d[~d.gene.isin(genes+['PTPN22'])].reset_index(drop=True);pool['pool_index']=np.arange(len(pool));save(pool,Q/f'REFERENCE_POOL_{system}.tsv');pools.append(pool)
 rng=np.random.default_rng(2026090831+i);indices=np.stack([rng.choice(len(pool),35,replace=False) for _ in range(10000)])
 save(pd.DataFrame(indices,columns=[f'g{j}' for j in range(35)]),Q/f'RANDOM_SET_INDICES_{system}.tsv')
 rr=pool.rank_score.to_numpy()[indices].mean(axis=1);refall.append(pd.DataFrame(dict(system=system,draw=np.arange(10000),program_shift=rr,seed=2026090831+i)))
 refsummary.append(reference_summary(rr,obs,system))
 save(d,Q/f'ALL_GENE_RANKS_{system}.tsv')
save(pd.concat(allranks,ignore_index=True),OUT/(P+'_GENE_LEVEL_ORIENTED_RANKS.tsv'));save(program,OUT/(P+'_SYSTEM_PROGRAM_SHIFT.tsv'));save(pd.concat(refall,ignore_index=True),OUT/(P+'_SYSTEM_RANDOM_REFERENCE.tsv'));save(refsummary,OUT/(P+'_SYSTEM_RANDOM_REFERENCE_SUMMARY.tsv'))
common=sorted(set.intersection(*[set(bg[bg.system.eq(s)].gene) for s in systems]));assert set(genes).issubset(common)
cd=pd.DataFrame(dict(gene=common));commonlong=[]
for system in systems:
 d=bg[bg.system.eq(system)].set_index('gene').loc[common];r=rank(d.oriented_effect.to_numpy());cd[system]=r
 commonlong.append(pd.DataFrame(dict(system=system,gene=common,native_effect=d.native_effect.to_numpy(),oriented_effect=d.oriented_effect.to_numpy(),common_rank_score=r)))
cd['CONSENSUS_PERTURBATION_SCORE']=cd[systems].mean(axis=1);cd['is_AF35']=cd.gene.isin(genes);cd['is_PTPN22']=cd.gene.eq('PTPN22')
save(pd.concat(commonlong,ignore_index=True),OUT/(P+'_COMMON_BACKGROUND.tsv'));save(cd,OUT/(P+'_CROSS_SYSTEM_GENE_CONSENSUS.tsv'))
observed=float(cd.set_index('gene').loc[genes,'CONSENSUS_PERTURBATION_SCORE'].mean());pool=cd[~cd.gene.isin(genes+['PTPN22'])].reset_index(drop=True);pool['pool_index']=np.arange(len(pool));assert len(pool)>=1000
save(pool,Q/'CROSS_SYSTEM_REFERENCE_POOL.tsv');rng=np.random.default_rng(2026090834);indices=np.stack([rng.choice(len(pool),35,replace=False) for _ in range(10000)]);save(pd.DataFrame(indices,columns=[f'g{j}' for j in range(35)]),Q/'CROSS_SYSTEM_RANDOM_SET_INDICES.tsv')
random=pd.DataFrame(dict(draw=np.arange(10000)))
for system in systems:random[system]=pool[system].to_numpy()[indices].mean(axis=1)
random['CROSS_SYSTEM_PROJECTION']=random[systems].mean(axis=1);random['seed']=2026090834
assert np.max(abs(random.CROSS_SYSTEM_PROJECTION.to_numpy()-pool.CONSENSUS_PERTURBATION_SCORE.to_numpy()[indices].mean(axis=1)))<1e-12
save(random,OUT/(P+'_CROSS_SYSTEM_RANDOM_REFERENCE.tsv'))
combined=reference_summary(random.CROSS_SYSTEM_PROJECTION.to_numpy(),observed,'COMMON_UNIVERSE_EQUAL_CONDITION_WEIGHT');combined['AF35_CROSS_SYSTEM_PROJECTION']=observed;combined['common_gene_n']=len(common);combined['reference_pool_n']=len(pool);combined['primary_condition_n']=3;combined['independent_experiment_n']=1;combined['seed']=2026090834
positive=sum(r['AF35_PROGRAM_SHIFT']>0 for r in program);strong_disagreement=any(r['observed']>r['random_q975'] for r in refsummary) and any(r['observed']<r['random_q025'] for r in refsummary)
state='COORDINATED_PROGRAM_SHIFT_SUPPORTED' if observed>0 and observed>combined['random_q975'] and positive>=2 else 'DIRECTIONALLY_SUPPORTIVE_NOT_EXCEPTIONAL' if observed>0 and observed<=combined['random_q975'] else 'PERTURBATION_CONTEXT_HETEROGENEOUS' if strong_disagreement else 'PROGRAM_LEVEL_SHIFT_NOT_SUPPORTED'
combined['scientific_state']=state;combined['secondary_concordance']='STRONG_CROSS_SYSTEM_DIRECTIONAL_CONCORDANCE' if positive==3 else 'NOT_ALL_PRIMARY_CONDITIONS_POSITIVE';save([combined],OUT/(P+'_CROSS_SYSTEM_PROJECTION_SUMMARY.tsv'))
cons=[]
for j,gene in enumerate(genes):
 row=dict(gene_rank=j+1,gene=gene);v=[]
 for d in allranks:
  r=d.set_index('gene').loc[gene];row[r['system']+'_oriented_effect']=r.oriented_effect;row[r['system']+'_rank_score']=r.rank_score;v.append(r.oriented_effect)
 row['positive_system_n']=int(np.sum(np.array(v)>0));row['negative_system_n']=int(np.sum(np.array(v)<0));row['zero_system_n']=int(np.sum(np.array(v)==0));row['positive_all']=row['positive_system_n']==3;row['positive_at_least_two']=row['positive_system_n']>=2;row['consistently_negative']=row['negative_system_n']==3;row['discordant']=row['positive_system_n']>0 and row['negative_system_n']>0;cons.append(row)
save(cons,OUT/(P+'_DIRECTIONAL_CONSISTENCY.tsv'))
js(Q/'PROJECTION_COMPLETION.json',dict(utc=now(),scientific_state=state,clinical_data_accessed=False,positive_condition_n=positive,strong_reference_disagreement=strong_disagreement,random_sets_per_reference=10000,reference_n=4))
save(access,Q/'PROJECTION_FILE_ACCESS_LOG.tsv')
print('Projection artifacts completed; inspect saved summaries after independent verification.',flush=True)
