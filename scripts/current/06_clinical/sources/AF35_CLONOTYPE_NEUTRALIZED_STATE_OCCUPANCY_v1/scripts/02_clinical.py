from common import *
from itertools import combinations
from math import comb
lock=jl(OUT/(P+'_FREEZE_LOCK.json'));assert not lock['outcome_loaded'] and not (OUT/(P+'_OUTCOME_JOIN_TIMESTAMP.txt')).exists()
for n,h in lock['hashes'].items():assert sha(OUT/n)==h,n
for n,h in lock['source_hashes'].items():assert sha(n)==h,n
stamp=now();write(P+'_OUTCOME_JOIN_TIMESTAMP.txt',stamp);js('QA/OUTCOME_OPEN_EVENT.json',dict(utc=stamp,freeze_utc=lock['utc'],freeze_lock_sha256=sha(OUT/(P+'_FREEZE_LOCK.json')),all_frozen_hashes_verified=True))
assert stamp>lock['utc']
label_source=K/'inputs/RESPONSE_LABEL_MAP_AFTER_FREEZE.tsv';km=read(K/'COMPLETE_SHA256_MANIFEST.tsv');kh=km.set_index('path').loc['inputs\\RESPONSE_LABEL_MAP_AFTER_FREEZE.tsv','sha256'];assert sha(label_source)==kh
labels=read(label_source,usecols=['patient','response_group','response']);assert labels.patient.is_unique and set(labels.response_group)=={'R','NR'} and (labels.response==labels.response_group.map({'R':1,'NR':0})).all()
source_labels=read(D/'authority/03_processed_data/phase3/GSE235863_BARCODE_CLONOTYPE_STATE_PTPN22.tsv.gz',usecols=['patient','response_group']).drop_duplicates();assert source_labels.patient.is_unique
assert labels.set_index('patient').response_group.equals(source_labels.set_index('patient').loc[labels.patient,'response_group'])
save(labels,'inputs/RESPONSE_LABEL_MAP_AFTER_FREEZE.tsv')
pc=read(OUT/(P+'_PATIENT_COMPONENTS_RESPONSE_BLIND.tsv'));cov=read(OUT/(P+'_TCR_COVERAGE.tsv'));d=pc.merge(cov[['patient','molecular_patient_eligible']],on='patient',validate='one_to_one').merge(labels,on='patient',validate='one_to_one',how='left');d['clinical_eligible']=d.molecular_patient_eligible&d.response.notna();save(d,'inputs/PATIENT_COMPONENTS_AFTER_RESPONSE_JOIN.tsv')
eligible=d[d.clinical_eligible].sort_values('patient');n=len(eligible);r=eligible.response.to_numpy()==1;nr=int(r.sum());nn=int((~r).sum());passed=n>=6 and min(nr,nn)>=2
save([dict(K_patients=len(cov),evaluable_patients=n,R_n=nr,NR_n=nn,resolved_cells=int(eligible.resolved_cells.sum()),clonotypes=int(eligible.clone_count.sum()),status='PASS' if passed else 'NOT_EVALUABLE_LOW_TCR_COVERAGE')],P+'_FEASIBILITY.tsv')
assert passed,'NOT_EVALUABLE_LOW_TCR_COVERAGE; stop primary inference'
patients=eligible.patient.tolist();keys=['O_CELL','O_CLONE','D_SIZE'];x=eligible[keys].to_numpy();delta=x[r].mean(axis=0)-x[~r].mean(axis=0);balance=delta[1]-delta[2]
def classify(v):
 cell,clone,size=v
 if cell<=0:return 'TCR_SUBSET_DOES_NOT_REPRODUCE_POSITIVE_OCCUPANCY'
 if clone>0:return 'MIXED_REPERTOIRE_AND_CLONE_SIZE' if size>0 else 'CLONE_NEUTRALIZED_OCCUPANCY_RETAINED'
 assert size>0
 return 'CLONE_SIZE_WEIGHTING_REQUIRED_FOR_POSITIVE_OCCUPANCY'
rng=np.random.default_rng(lock['bootstrap_seed']);ri=np.flatnonzero(r);ni=np.flatnonzero(~r);idx=np.column_stack([rng.choice(ri,size=(10000,nr),replace=True),rng.choice(ni,size=(10000,nn),replace=True)])
save(pd.DataFrame(idx,columns=[f'slot_{i+1}' for i in range(n)]),'inputs/BOOTSTRAP_INDICES_ZERO_BASED.tsv');save({'position_zero_based':range(n),'patient':patients,'response':r.astype(int)},'inputs/BOOTSTRAP_PATIENT_ORDER.tsv')
bd=x[idx[:,:nr]].mean(axis=1)-x[idx[:,nr:]].mean(axis=1);bb=bd[:,1]-bd[:,2];boot=pd.DataFrame(bd,columns=['Delta_'+k for k in keys]);boot['Delta_BALANCE']=bb;boot.insert(0,'draw',range(1,10001));save(boot,P+'_BOOTSTRAP.tsv')
effects=[];summary=[]
for j,k in enumerate(keys+['BALANCE']):
 vv=bd[:,j] if j<3 else bb;obs=delta[j] if j<3 else balance;lo,hi=np.quantile(vv,[.025,.975])
 effects.append(dict(component=k,R_n=nr,NR_n=nn,mean_R=float(x[r,j].mean()) if j<3 else float((x[r,1]-x[r,2]).mean()),mean_NR=float(x[~r,j].mean()) if j<3 else float((x[~r,1]-x[~r,2]).mean()),raw_R_minus_NR=obs,bootstrap_ci_low=lo,bootstrap_ci_high=hi,units='FROZEN_K_ALPHA_AF35_UNITS'))
 summary.append(dict(component=k,draws=10000,observed=obs,bootstrap_median=float(np.median(vv)),ci_low=lo,ci_high=hi,interval='PATIENT_RESPONSE_STRATIFIED_PERCENTILE_95',frozen_components=True))
save(effects,P+'_CLINICAL_DECOMPOSITION.tsv');save(summary,P+'_BOOTSTRAP_SUMMARY.tsv')
lo,hi=np.quantile(bb,[.025,.975]);leading='REPERTOIRE_COMPONENT_LEADING' if lo>0 else ('CLONE_SIZE_COMPONENT_LEADING' if hi<0 else 'MIXED_OR_IMPRECISE');pattern=classify(delta)
kp=read(K/'AF35_STATE_OCCUPANCY_PATIENT_COMPONENTS_RESPONSE_BLIND.tsv').query("compartment == 'P'").merge(labels,on='patient',validate='one_to_one').sort_values('patient');kall=kp.loc[kp.response==1,'O_i'].mean()-kp.loc[kp.response==0,'O_i'].mean();ke=kp.set_index('patient').loc[patients,'O_i'].to_numpy();kdelta=ke[r].mean()-ke[~r].mean();bridgegate='POSITIVE_DIRECTION' if kdelta>0 and delta[0]>0 else 'DIRECTION_NOT_REPLICATED_IN_TCR_SUBSET'
bridge=read(OUT/(P+'_K_BRIDGE_RESPONSE_BLIND.tsv')).merge(labels,on='patient',validate='one_to_one');bridge['clinical_eligible']=bridge.patient.isin(patients);bridge['Delta_O_K_full_K']=kall;bridge['Delta_O_K_same_L_patients']=kdelta;bridge['Delta_O_CELL']=delta[0];bridge['universe_bridge_gate']=bridgegate;save(bridge,P+'_K_BRIDGE.tsv')
perms=[];assignments=[];pdeltas=[]
for j,cmb in enumerate(combinations(range(n),nr)):
 rr=np.zeros(n,dtype=bool);rr[list(cmb)]=True;pdv=x[rr].mean(axis=0)-x[~rr].mean(axis=0);pdeltas.append(pdv);assignments.append(dict(assignment_id=j+1,responder_patients=';'.join(patients[i] for i in cmb),Delta_O_CELL=pdv[0],Delta_O_CLONE=pdv[1],Delta_D_SIZE=pdv[2]))
pdeltas=np.array(pdeltas);assert len(pdeltas)==comb(n,nr)
for j,k in enumerate(keys):
 extreme=int((abs(pdeltas[:,j])>=abs(delta[j])-1e-12).sum());perms.append(dict(component=k,observed=delta[j],patient_n=n,R_n=nr,NR_n=nn,unique_assignments=len(pdeltas),extreme_assignments=extreme,exact_two_sided_P=extreme/len(pdeltas),supportive_only=True))
save(perms,P+'_EXACT_PERMUTATION.tsv');save(assignments,'QA/EXACT_PERMUTATION_ASSIGNMENTS.tsv')
prop=read(OUT/'AF35_CLONE_NEUTRALIZATION_PATIENT_STATE_PROPORTIONS_RESPONSE_BLIND.tsv');order=read(OUT/'inputs/RELEASED_STATE_ORDER.tsv');states=order.released_state.tolist();state_rows=[]
for state_order,s in enumerate(states,1):
 z=prop[prop.released_state==s].set_index('patient').loc[patients];a=float(z.alpha_s.iloc[0]);cv=(z.p_cell.to_numpy()[r].mean()-z.p_cell.to_numpy()[~r].mean())*a;ev=(z.p_clone.to_numpy()[r].mean()-z.p_clone.to_numpy()[~r].mean())*a
 state_rows.append(dict(state_order=state_order,released_state=s,alpha_s=a,mean_p_cell_R=z.p_cell.to_numpy()[r].mean(),mean_p_cell_NR=z.p_cell.to_numpy()[~r].mean(),mean_p_clone_R=z.p_clone.to_numpy()[r].mean(),mean_p_clone_NR=z.p_clone.to_numpy()[~r].mean(),CELL_WEIGHTED_STATE_CONTRIBUTION=cv,CLONE_EQUAL_STATE_CONTRIBUTION=ev,CLONE_SIZE_STATE_CONTRIBUTION=cv-ev))
st=pd.DataFrame(state_rows);save(st,P+'_STATE_CONTRIBUTIONS.tsv')
state_error=st[['CELL_WEIGHTED_STATE_CONTRIBUTION','CLONE_EQUAL_STATE_CONTRIBUTION','CLONE_SIZE_STATE_CONTRIBUTION']].sum().to_numpy()-delta;err=float(abs(delta[0]-delta[1:].sum()));berr=float(np.max(abs(bd[:,0]-bd[:,1:].sum(axis=1))));assert max(err,berr,np.max(abs(state_error)))<1e-12
save([dict(group_identity_error=err,bootstrap_identity_max_error=berr,state_CELL_sum_error=abs(state_error[0]),state_CLONE_sum_error=abs(state_error[1]),state_SIZE_sum_error=abs(state_error[2]),status='PASS')],'QA/CLINICAL_EXACT_IDENTITY_VERIFICATION.tsv')
pair=read(OUT/'inputs/PAIRED_CDR3_PATIENT_COMPONENTS_RESPONSE_BLIND.tsv').set_index('patient').loc[patients];px=pair[keys].to_numpy();pdelta=px[r].mean(axis=0)-px[~r].mean(axis=0);paired_pattern=classify(pdelta)
sens=[]
for i,pat in enumerate(patients):sens.append(dict(row_type='PATIENT',patient=pat,analysis='PAIRED_CDR3_IDENTITY_SENSITIVITY',author_clones=int(eligible.iloc[i].clone_count),paired_clones=int(pair.iloc[i].clone_count),author_O_CELL=x[i,0],author_O_CLONE=x[i,1],author_D_SIZE=x[i,2],paired_O_CELL=px[i,0],paired_O_CLONE=px[i,1],paired_D_SIZE=px[i,2],O_CLONE_change=px[i,1]-x[i,1],D_SIZE_change=px[i,2]-x[i,2]))
sens.append(dict(row_type='R_MINUS_NR',patient='AGGREGATE',analysis='PAIRED_CDR3_IDENTITY_SENSITIVITY',author_O_CELL=delta[0],author_O_CLONE=delta[1],author_D_SIZE=delta[2],paired_O_CELL=pdelta[0],paired_O_CLONE=pdelta[1],paired_D_SIZE=pdelta[2],O_CLONE_change=pdelta[1]-delta[1],D_SIZE_change=pdelta[2]-delta[2],scientific_pattern=paired_pattern));save(sens,P+'_PAIRED_CDR3_SENSITIVITY.tsv')
dec=dict(scientific_pattern=pattern,balance_pattern=leading,Delta_BALANCE=balance,balance_ci_low=lo,balance_ci_high=hi,TCR_universe_bridge=bridgegate,paired_identity_pattern=paired_pattern)
save([dec],P+'_SCIENTIFIC_INTERPRETATION.tsv')
write(P+'_OUTCOME_FIREWALL_AUDIT.md',f'''# L 结局防火墙

全部患者/克隆分量、K bridge、状态比例、耦合、集中度及配对CDR3敏感性分量于{lock['utc']}冻结。首次响应列读取事件为{stamp}，严格在后；事件JSON记录当时冻结锁SHA256。接入前验证全部{len(lock['hashes'])}个冻结文件与{len(lock['source_hashes'])}个直接源文件哈希。
响应来自K原映射，另与Clone-State原response_group逐患者核对；本轮无响应参与TCR解析、local clone size、f、A、O或D构造。K已有结果及用户背景不等于本轮绝对未知，但本轮执行隔离、固定公式和时序可核验。仅在冻结后判断响应覆盖和组人数，并运行同一患者集合的比较。
所有local clone counts只来自K pre P精确细胞集合，未读global clone.size、post/P-T联合大小、持续/后续观察类别。K冻结alpha原样复用，无重新拟合。既有配对CDR3敏感性同样响应前冻结，不替换主身份。
''')
js('QA/CLINICAL_COMPLETION_LOCK.json',dict(utc=now(),patient_n=n,R_n=nr,NR_n=nn,scientific_decision=dec,hashes={p.name:sha(p) for p in OUT.glob(P+'_*.tsv')}))
print(pd.DataFrame(effects).to_string(index=False));print(pd.DataFrame([dec]).to_string(index=False));print(pd.DataFrame(perms).to_string(index=False));print('Paired contrasts',pdelta)
