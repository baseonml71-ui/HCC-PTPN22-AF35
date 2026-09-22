from common import *
from matplotlib.patches import Rectangle
from PIL import Image
import xml.etree.ElementTree as ET
import base64,io
q=Work(4,600,730)
q.placeholder('A',8,20,544,26)
palette=q.src('F3_STATE_COLOR_DICTIONARY_B3.tsv');states=[r['state_id'] for r in sorted(palette,key=lambda r:int(r['display_order']))];colors={r['state_id']:r['hex_color'] for r in palette}
scores=q.src('AF35_STATE_FROZEN_SCORE_SUMMARY.tsv');program=q.src('GSE235863_STATE_PROGRAM_MATRIX.tsv');comp=q.src('GSE235863_PATIENT_CLONE_CONTEXT_STATE_COMPOSITION.tsv')
dp=q.src('F4D_THREE_STAGES_PATIENTS.tsv');ds=q.src('CLONE_STATE_EFFECT_SUMMARY.tsv');seven=q.src('F4E_SEVEN_FROZEN_COMPONENTS.tsv');paired=q.src('AF35_CLONE_NEUTRALIZATION_PATIENT_COMPONENTS_RESPONSE_BLIND.tsv');corr=q.src('AF35_CLONE_SIZE_STATE_COUPLING.tsv');times=q.src('AF35_M2_PATIENT_TIMEPOINT_COMPONENTS_RESPONSE_BLIND.tsv');theta=q.src('AF35_M2_CLINICAL_LONGITUDINAL_CONTRASTS.tsv');overlap=q.src('CLONE_CONTEXT_ELIGIBILITY_AWARE_OVERLAP.tsv');triple=q.src('CLONE_CONTEXT_SHARED_UNIVERSE_TRIPLE_OVERLAP.tsv')
assert len(states)==20 and len(dp)==23 and len(paired)==7 and len(times)==14
q.panel('A',12,14,276,191,'Frozen T-cell atlas')
atlas=ROOT/'绘图视觉/HCC_ICI_FIGURE3_V2_1B3_ATLAS_FINALIZATION_v1'
svg=atlas/'F3B_B3_1_DIRECT_LABEL_ONLY.svg';node=ET.parse(svg).find('.//{http://www.w3.org/2000/svg}image');blob=base64.b64decode(node.attrib['{http://www.w3.org/1999/xlink}href'].split(',')[1]);im=Image.open(io.BytesIO(blob))
ax=q.ax(22,35,255,166);ax.imshow(im,extent=[-10.788582293578,15.816922293578,-5.24886850152905,12.06],interpolation='none');ax.set_aspect('equal');ax.axis('off')
labels=read(atlas/'F3B_B3_LABEL_PLACEMENT.tsv')
atlas_manifest=json.loads((atlas/'F3B_B3_MANIFEST.json').read_text(encoding='utf-8-sig'))
for name in ['F3B_B3_1_DIRECT_LABEL_ONLY.svg','F3B_B3_LABEL_PLACEMENT.tsv']:
    assert sha(atlas/name)==next(r['sha256'] for r in atlas_manifest['files'] if r['relative_path']==name)
for r in labels:
    x,y=float(r['display_label_x']),float(r['display_label_y']);cx,cy=float(r['source_x_center']),float(r['source_y_center'])
    ax.annotate(r['state_label'],(cx,cy),xytext=(x,y),ha='center',va='center',fontsize=11,bbox=dict(facecolor='white',edgecolor='none',pad=.3),arrowprops=dict(arrowstyle='-',color=GRAY,lw=.6) if r['leader_line_used']=='YES' else None)
    q.mark('A',state=r['state_id'],label_x=x,label_y=y,color=colors[r['state_id']])
for p in [svg,atlas/'F3B_B3_LABEL_PLACEMENT.tsv']:q.binds.append(dict(panel='A frozen visual',path=str(p),sha256=sha(p)))
q.mark('A',kind='frozen_scatter_raster',sha256=hashlib.sha256(blob).hexdigest(),native_pixels=list(im.size),unchanged_embedding=True)

q.panel('B',310,14,278,191,'State AF35 and functional localization')
progs=list(dict.fromkeys(r['program'] for r in program));assert len(progs)==12
ax=q.ax(371,44,212,150)
arr=np.array([[float(next(r for r in program if r['released_t_cell_state']==s and r['program']==p)['patient_mean_score']) for p in progs] for s in states])
lim=np.ceil(np.max(np.abs(arr))*10)/10
im=ax.pcolormesh(np.arange(13)-.5,np.arange(21)-.5,arr,cmap=DIV,vmin=-lim,vmax=lim);ax.invert_yaxis();ax.set_xticks(range(12),progs,rotation=90);ax.xaxis.tick_top();ax.set_yticks(range(20),states);ax.tick_params(length=0);ax.spines[['left','bottom']].set_visible(False)
for tick in ax.get_yticklabels():tick.set_color(INK)
ax.tick_params(axis='y',pad=10)
from matplotlib.colors import ListedColormap
stripe=q.ax(367,44,2,150);stripe.pcolormesh([0,1],np.arange(21),np.arange(20).reshape(-1,1),cmap=ListedColormap([colors[s] for s in states]));stripe.invert_yaxis();stripe.axis('off')
for i,s in enumerate(states):
    for j,p in enumerate(progs):
        q.mark('B',state=s,program=p,value=arr[i,j])
# Frozen AF35 kept on a separate native numeric rail, never on program z axis.
for i,s in enumerate(states):
    r=next(r for r in scores if r['scope']=='ALL_MAPPED' and r['released_state']==s)
    q.text(310,47+i*7.5,f"{float(r['mean_frozen_AF35']):.3f}  {r['patient_n']}",9);q.mark('B',kind='AF35',**r)
q.text(310,35,'AF35  n',10,True)
q.fig.colorbar(im,cax=q.ax(390,201,155,2),orientation='horizontal',ticks=[-lim,0,lim]).set_label('Program score',fontsize=9)

q.panel('C',12,230,276,230,'Clone-context state occupancy')
contexts=list(dict.fromkeys((r['context_family'],r['context_level']) for r in comp));assert len(contexts)==6
frozen_composition_path=ROOT/'绘图视觉/HCC_ICI_FIGURE3_V2_1C_HEATMAP_CALIBRATION_v1/qa/F3D_EQUAL_PATIENT_DISPLAY_VALUES.tsv'
frozen_composition=read(frozen_composition_path)
composition_manifest=json.loads((frozen_composition_path.parent.parent/'F3D_F3F_MANIFEST.json').read_text(encoding='utf-8-sig'))
assert sha(frozen_composition_path)==next(r['sha256'] for r in composition_manifest['files'] if r['relative_path']=='qa/F3D_EQUAL_PATIENT_DISPLAY_VALUES.tsv')
q.binds.append(dict(panel='C existing equal-patient display',path=str(frozen_composition_path),sha256=sha(frozen_composition_path)))
shutil.copy2(frozen_composition_path,q.ready/frozen_composition_path.name)
ax=q.ax(70,265,213,165);arr=np.zeros((20,6));ns=[]
for j,(fam,lev) in enumerate(contexts):
    rr=[r for r in comp if r['context_family']==fam and r['context_level']==lev];ns.append(len(set(r['patient'] for r in rr)))
    for i,s in enumerate(states):
        vals=[float(r['state_proportion_within_patient_context']) for r in rr if r['sub_cluster']==s]
        frozen=next(r for r in frozen_composition if (r['context_family'],r['context_level'],r['sub_cluster'])==(fam,lev,s))
        arr[i,j]=float(frozen['value']);assert abs(arr[i,j]-np.mean(vals))<1e-12
        q.mark('C',family=fam,level=lev,state=s,patient_n=len(vals),display_patient_mean=arr[i,j],source_values=vals)
assert ns==[9,9,7,7,7,7]
ax.pcolormesh(np.arange(7)-.5,np.arange(21)-.5,arr,cmap=SEQ,vmin=0,vmax=max(.3,arr.max()));ax.invert_yaxis();ax.set_yticks(range(20),states);ax.set_xticks(range(6),[r[1] for r in contexts],rotation=55,ha='left');ax.xaxis.tick_top();ax.tick_params(length=0);ax.spines[['left','bottom']].set_visible(False)
for tick in ax.get_yticklabels():tick.set_color(INK)
ax.tick_params(axis='y',pad=10)
stripe=q.ax(66,265,2,165);stripe.pcolormesh([0,1],np.arange(21),np.arange(20).reshape(-1,1),cmap=ListedColormap([colors[s] for s in states]));stripe.invert_yaxis();stripe.axis('off')
for i in range(20):
    for j in range(6):ax.text(j,i,f'{arr[i,j]*100:.1f}',ha='center',va='center',fontsize=9,color='white' if arr[i,j]>.18 else INK)
q.text(70,437,'Patient-mean state composition (%)',10)
q.text(70,443,'Expansion 9 / Persistence 7 / Later 7 patients',8)
# Mathematical relation applies only to the explicitly named eligible sets.
q.text(25,449,r'Eligible: Persistence $\subseteq$ Expanded; Later $\subseteq$ Expanded',8)
q.text(25,456,'Common baseline-PB positive sets: Persistence = Later',8)
assert all(int(r['n01'])==0 for r in overlap if r['status']=='JOINTLY_ELIGIBLE' and r['context1']=='expanded')
q.mark('C',kind='structural_dependence',records=overlap,shared_baseline_pb=triple)

q.panel('D',310,230,278,230,'Patient analytical-condition contrasts')
fields=['raw_effect','supported_cell_raw_effect','adjusted_effect'];labels=['Original','Same-support','State-adjusted']
for ci,(context,label,n) in enumerate([('EXPANSION','Expansion',9),('PERSISTENCE','Persistence',7),('LATER_OBSERVATION','Later-observed',7)]):
    rs=[r for r in dp if r['context']==context];assert len(rs)==n
    top=[259,332,393][ci];ax=q.ax(370,top,205,n*5.7);a=np.array([[float(r[f]) for f in fields] for r in rs]);assert np.abs(a).max()<.4
    ax.pcolormesh(np.arange(4)-.5,np.arange(n+1)-.5,a,cmap=DIV,vmin=-.4,vmax=.4);ax.invert_yaxis();ax.set_yticks(range(n),[r['patient'] for r in rs]);ax.set_xticks(range(3),labels if ci==0 else ['']*3);ax.xaxis.tick_top();ax.tick_params(length=0);ax.spines[['left','bottom']].set_visible(False)
    q.text(310,top+4,label,10,True)
    for i,r in enumerate(rs):
        for j,f in enumerate(fields):ax.text(j,i,f'{a[i,j]:+.3f}',ha='center',va='center',fontsize=10);q.mark('D',context=context,patient=r['patient'],field=f,value=a[i,j])
    for j,est in [(0,'RAW'),(2,'STATE_ADJUSTED')]:
        r=next(r for r in ds if r['identity']=='AUTHOR_CLONE_ID' and r['score']=='AF35' and r['context']==context and r['estimate']==est)
        q.text(372+j*68,top+n*5.7+5,f"{float(r['mean_patient_effect']):+.3f} [{float(r['bootstrap_ci_low']):+.3f}, {float(r['bootstrap_ci_high']):+.3f}]",9);q.mark('D',kind='summary',**r)

q.panel('E',12,485,300,224,'Signed K and L decomposition')
keys=['M','O','O_CLONE','D_SIZE','W','X','O_minus_W'];labs=['Overall','Occupancy','Clone-equal','Clone-size','Within-state','Interaction','ΔO − ΔW'];ax=q.ax(88,506,208,180)
for i,(k,label) in enumerate(zip(keys,labs)):
    r=next(r for r in seven if r['component']==k);v,a,b=map(float,[r['raw_R_minus_NR'],r['bootstrap_ci_low'],r['bootstrap_ci_high']]);y=6-i
    ax.barh(y,b-a,left=a,height=.52,color=GRAY,alpha=.17);ax.barh(y,v,height=.19,color=TEAL if v>=0 else BLUE);ax.text(.203,y,f'{v:+.3f}',va='center',fontsize=10);q.mark('E',**r)
ax.axvline(0,color=GRAY,lw=.65);ax.set_yticks(range(6,-1,-1),labs);ax.set(xlim=(-.14,.25),ylim=(-.6,6.6),xlabel='AF35 difference (R − NR)');ax.spines[['left','bottom']].set_visible(False)

q.panel('F',330,485,258,77,'Clonotype weighting within patients')
ax=q.ax(369,503,191,48)
for i,r in enumerate(paired):
    a,b=map(float,[r['O_CELL'],r['O_CLONE']]);ax.plot([a,b],[6-i]*2,color=GRAY,lw=.8);ax.scatter(a,6-i,s=20,c=TEAL);ax.scatter(b,6-i,s=20,facecolor='white',edgecolor=TEAL)
    cr=next(z for z in corr if z['patient']==r['patient']);q.mark('F',**r);q.mark('F',kind='stored_coupling',**cr)
    q.text(563,507+i*6.86,f"{float(cr['spearman_log1p_size_potential']):+.3f}",8)
ax.set_yticks(range(6,-1,-1),[r['patient'] for r in paired]);ax.set_xlabel('AF35 representation');ax.spines[['left','bottom']].set_visible(False);q.text(369,496,'● Cell-weighted     ○ Clone-equal',10)
q.text(563,496,'ρ',9)

q.panel('G',330,584,258,68,'Paired blood representation')
for j,f in enumerate(['O_CELL','O_CLONE','D_SIZE']):
    ax=q.ax(356+j*78,600,63,43)
    endpoints=[]
    for i,p in enumerate([r['patient'] for r in paired]):
        rs=[next(r for r in times if r['patient']==p and r['timepoint']==t) for t in ['pre','post']];vs=[float(r[f]) for r in rs];ax.plot([0,1],vs,color=GRAY,lw=.6,marker='o',ms=2);endpoints.append((vs[1],p))
        for r in rs:q.mark('G',field=f,**r)
    lower,upper=ax.get_ylim();step=(upper-lower)*.095;last=lower
    for value,p in sorted(endpoints):
        label_y=max(value,last+step);last=label_y
        ax.annotate(p,(1,value),xytext=(1.12,label_y),fontsize=9,va='center',arrowprops=dict(arrowstyle='-',lw=.4,color=GRAY))
    ax.set_ylim(lower,max(upper,last+step));ax.set_xticks([0,1],['Pre','Later']);ax.set_xlim(-.1,1.5);ax.set_title(['Cell-weighted','Clone-equal','Clone-size'][j],fontsize=10);ax.spines[['left','bottom']].set_visible(False)

q.panel('H',330,681,258,40,'Global longitudinal contrast')
ax=q.ax(387,686,160,29)
for i,r in enumerate(theta[:3]):
    v,a,b=map(float,[r['estimate'],r['bootstrap_ci_low'],r['bootstrap_ci_high']]);ax.barh(2-i,b-a,left=a,height=.45,color=GRAY,alpha=.18);ax.barh(2-i,v,height=.14,color=BLUE);q.text(552,689+i*9,f"P {float(r['exact_two_sided_p']):.3f}",9);q.mark('H',**r)
ax.axvline(0,color=GRAY,lw=.65);ax.set_yticks([2,1,0],['θ cell','θ clone','θ size']);ax.set(xlim=(-.08,.06),ylim=(-.5,2.5));ax.set_xticks([-.05,0,.05]);ax.spines[['left','bottom']].set_visible(False)
q.save()
