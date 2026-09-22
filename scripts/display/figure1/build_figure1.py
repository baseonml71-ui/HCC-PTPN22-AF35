"""Deterministic display of frozen Figure 1 evidence. No statistical estimation."""
from pathlib import Path
import csv, json, hashlib, shutil, math, importlib.util
import pymupdf as fitz
import numpy as np
from PIL import Image

ROOT=Path('.')
OUT=Path(__file__).resolve().parent
PAY=ROOT/'HCC_ICI_MAIN_FIGURE_FROZEN_SOURCE_PAYLOAD_V2'
QA=OUT/'qa'; READY=OUT/'figure_ready'; COPIES=OUT/'source_tables'
for folder in [QA,READY,COPIES]:folder.mkdir(exist_ok=True)
MM=72/25.4
spec=importlib.util.spec_from_file_location('drawing_only',ROOT/'绘图视觉/HCC_ICI_NEW_HIGH_VALUE_RESULT_PANELS_v1/qa/build_panels.py')
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base);base.QA=QA
D=base.Drawing
INK,GRAY,LIGHT=base.INK,base.GRAY,base.LIGHT
BLUE,RED,TEAL=base.BLUE,base.RED,base.TEAL
WHITE='#FFFFFF'
PREFIX='FIGURE1_FINAL_DATA_DRIVEN_V1'

def read(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def table(p,rows):
    with Path(p).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def dump(p,data):Path(p).write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
def num(r,k):
    try:return float(r[k])
    except (ValueError,TypeError):return math.nan
def signed(v,digits=3):return f'{v:+.{digits}f}' if math.isfinite(v) else 'NA'

bindings=[]; sources={}
def bind(panel,p,h):
    p=Path(p);assert sha(p)==h,(p,'FROZEN_HASH_MISMATCH')
    dest=COPIES/p.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True)
    if not dest.exists():shutil.copy2(p,dest)
    assert sha(dest)==h
    sources[p.name]=p
    bindings.append(dict(panel=panel,path=str(p),copy=dest.relative_to(OUT).as_posix(),sha256=h))
    return p

ledger=read(PAY/'MAIN_FIGURE_39_PANEL_SOURCE_BINDING_LEDGER.tsv')
for row in ledger:
    if row['figure']!='1':continue
    for p,h in zip(row['authoritative_source_file'].split(';'),row['source_file_sha256'].split(';')):bind(row['panel'],p,h)
payload_manifest={r['path']:r['sha256'] for r in read(PAY/'DELIVERABLE_MANIFEST_SHA256.tsv')}
for rel in ['Figure1/display_support/F1B_CURRENT_PTPN22_STATE_COMPARISON.tsv','Figure1/display_support/F1D_PRIMARY_META_ONLY.tsv']:
    bind('frozen display support',PAY/rel,payload_manifest[rel])

def extra(panel,folder,rel,manifest):
    mp={r['path'].replace('\\','/'):r['sha256'] for r in read(ROOT/folder/manifest)}
    return bind(panel,ROOT/folder/rel,mp[rel])
for rel in ['07_robustness/SOURCE_GATE_SUMMARY.tsv','07_robustness/PTPN22_AF35_TECHNICAL_SENSITIVITY.tsv','02_manifest/PTPN22_AF35_COUPLING_ANALYSIS_MANIFEST_V1.md']:
    extra('1E/1F','PTPN22_AF35_COUPLING_ARCHITECTURE_v1',rel,'09_qa/DELIVERABLE_MANIFEST_SHA256.tsv')
for rel in ['11_checkpoint/PTPN22_STATE_AXIS_SPECIFICITY_CHECKPOINT.md','02_manifest/PTPN22_STATE_AXIS_SPECIFICITY_ANALYSIS_MANIFEST_V1.md']:
    extra('1F','PTPN22_STATE_AXIS_SPECIFICITY_v1',rel,'09_qa/DELIVERABLE_MANIFEST_SHA256.tsv')
extra('1G','PTPN22_AF35_PROGRAM_PERTURBATION_PROJECTION_v1','PTPN22_AF35_GENE_LEVEL_ORIENTED_RANKS.tsv','COMPLETE_SHA256_MANIFEST.tsv')
for p in [ROOT/'绘图视觉/HCC_ICI_NEW_HIGH_VALUE_RESULT_PANELS_v1/qa/build_panels.py',Path('private_dependencies/unavailable_display_context.txt')]:
    bindings.append(dict(panel='implementation / authorization',path=str(p),copy='',sha256=sha(p)))
table(OUT/'SOURCE_FILES_SHA256.tsv',bindings)

baseline=QA/'PROTECTED_BASELINE.json'
if not baseline.exists():
    protected={str(PAY/rel):h for rel,h in payload_manifest.items()}
    protected.update({r['path']:r['sha256'] for r in bindings})
    for folder in ['绘图视觉/VISUAL_GRAMMAR_PRODUCTION_PROOF_V1','绘图视觉/HCC_ICI_FIVE_FIGURE_FINAL_ASSEMBLY_PHASE_A_V1','绘图视觉/HCC_ICI_FIGURE3_SCIENTIFIC_SYNC_V1R1']:
        protected.update({str(p):sha(p) for p in (ROOT/folder).rglob('*') if p.is_file()})
    protected.update({str(p):sha(p) for p in ROOT.glob('*MANUSCRIPT*.md')})
    dump(baseline,protected)

def rows(name):return read(sources[name])
def ready(name,data):table(READY/(name+'.tsv'),data);return data
A=ready('A_MARKER_STATE',rows('FIGURE1_CANONICAL_MARKER_STATE_PLOT_READY.tsv'))
B=ready('B_STATE_RECURRENCE',rows('F1B_CURRENT_PTPN22_STATE_COMPARISON.tsv'))
C=ready('C_PROGRAM_STATE',rows('FIGURE1_STATE_PROGRAM_ASSOCIATION_PLOT_READY.tsv'))
COHORT=ready('D_FOUR_COHORTS',rows('FIGURE2_STRICT_FOUR_COHORT_EFFECTS.tsv'))
META=ready('D_PRIMARY_META',rows('F1D_PRIMARY_META_ONLY.tsv'))[0]
PRIMARY=rows('PTPN22_AF35_PATIENT_COUPLING_DECOMPOSITION.tsv')
TTECH=[r for r in rows('PTPN22_AF35_TECHNICAL_SENSITIVITY.tsv') if r['analysis']=='TECHNICAL_ADJUSTED']
LTECH=rows('PTPN22_AF35_LINEAGE_TECHNICAL_ADJUSTMENT.tsv')
patient_ids=[r['patient'] for r in PRIMARY]
E=[]
for lineage,org,tech in [('All T cells',PRIMARY,TTECH),('CD4',rows('PTPN22_AF35_CD4_COUPLING_DECOMPOSITION.tsv'),[r for r in LTECH if r['lineage']=='CD4']),('CD8',rows('PTPN22_AF35_CD8_COUPLING_DECOMPOSITION.tsv'),[r for r in LTECH if r['lineage']=='CD8'])]:
    for condition,data in [('Original',org),('Technical-adjusted',tech)]:
        lookup={r['patient']:r for r in data};assert set(lookup)==set(patient_ids)
        for patient in patient_ids:
            for component in ['C_total','C_between','C_within']:
                E.append(dict(lineage=lineage,patient=patient,condition=condition,component=component,value=lookup[patient][component]))
ready('E_PATIENT_COMPONENTS',E)
FT=rows('SOURCE_GATE_SUMMARY.tsv');FL=rows('LINEAGE_GATE_SUMMARY.tsv')
F=[]
for lineage in ['All T cells','CD4','CD8']:
    for condition in ['PRIMARY','TECHNICAL_ADJUSTED']:
        r=next(r for r in (FT if lineage=='All T cells' else FL) if r['analysis']==condition and (lineage=='All T cells' or r['lineage']==lineage))
        F.append(dict(lineage=lineage,condition=condition,median_C_total=r['median_C_total'],median_C_between=r['median_C_between'],median_C_within=r['median_C_within'],n=r['eligible_patient_n'],final_gate=r['final_gate']))
ready('F_LINEAGE_SUMMARY',F)
LOSO=ready('F_SINGLE_STATE_LOSO',[r for r in rows('PTPN22_AF35_LOSO_SUMMARY.tsv') if r['removal_type']=='SINGLE_STATE'])
G1=ready('G1_ORIENTED_GENE_EFFECTS',rows('PTPN22_AF35_GENE_LEVEL_ORIENTED_RANKS.tsv'))
GSH=ready('G1_PROGRAM_SHIFTS',rows('PTPN22_AF35_SYSTEM_PROGRAM_SHIFT.tsv'))
GR=ready('G1_RANDOM_REFERENCES',rows('PTPN22_AF35_SYSTEM_RANDOM_REFERENCE_SUMMARY.tsv'))
G2=ready('G2_DONORS',rows('GSE94859_AF35_ORDERED_DONOR_DELTAS.tsv'))
GS=ready('G2_GENOTYPE_EFFECT',rows('GSE94859_AF35_GENOTYPE_INTERACTION.tsv'))[0]
assert (len(A),len(B),len(C),len(E),len(LOSO),len(G1),len(G2))==(192,6,72,108,20,105,16)
assert sum(r['genotype']=='CC' for r in G2)==sum(r['genotype']=='TT' for r in G2)==8
assert len(COHORT)==4 and META['model_id']=='PRIMARY_BULK_TUMOR_REML_HKSJ'

states=[r['state'] for r in sorted({r['state']:r for r in A}.values(),key=lambda r:int(r['state_order']))]
state_labels=['CD8 progenitor-exh.','CD8 effector','CD8 terminal-exh.','Proliferating T','CD4 Tfh-like','Treg']
markers=[r['marker'] for r in sorted({r['marker']:r for r in A}.values(),key=lambda r:int(r['marker_order']))]
programs=[r['program'] for r in sorted({r['program']:r for r in C}.values(),key=lambda r:int(r['program_order']))]
program_labels=['TCR activation','NFAT','AP1','NFκB','IFN response','IFNγ / TNF','Cytotoxicity','Exhaustion','Checkpoint','Proliferation','Effector differentiation','Treg']

def color(v,limit):return base.blend(v,limit)
def bar(d,x,y,w,v,lo,hi,h=3,col=BLUE):
    assert lo<=v<=hi
    zero=x+w*(0-lo)/(hi-lo);xx=x+w*(v-lo)/(hi-lo)
    d.rect(min(zero,xx),y-h/2,max(abs(xx-zero),.08),h,col)
def haxis(d,x,y,w,lo,hi,ticks,label=None):
    d.line(x,y,x+w,y,GRAY,.22)
    for v in ticks:
        xx=x+w*(v-lo)/(hi-lo);d.line(xx,y,xx,y+1.2,GRAY,.2);d.text(xx,y+5,f'{v:g}',10,anchor='middle',color=GRAY)
    if label:d.text(x+w/2,y+11,label,10,anchor='middle')
def cbar(d,x,y,w,limit,label):
    for i in range(100):d.rect(x+i*w/100,y,w/100+.01,2.4,color((2*i/99-1)*limit,limit))
    for v in [-limit,0,limit]:d.text(x+w*(v/limit+1)/2,y+6.5,f'{v:g}',10,anchor='middle')
    d.text(x,y+12,label,10)
def cell(d,x,y,w,h,v,lim,number=True,outline=False):
    if not math.isfinite(v):
        d.rect(x,y,w,h,'#E5E7E8');d.line(x+1,y+1,x+w-1,y+h-1,GRAY,.3);d.line(x+1,y+h-1,x+w-1,y+1,GRAY,.3)
        return
    d.rect(x,y,w-.3,h-.3,color(v,lim))
    if outline:
        for a,b,c,e in [(x+.5,y+.5,x+w-.8,y+.5),(x+.5,y+h-.8,x+w-.8,y+h-.8),(x+.5,y+.5,x+.5,y+h-.8),(x+w-.8,y+.5,x+w-.8,y+h-.8)]:d.line(a,b,c,e,INK,.28)
    if number:d.text(x+w/2,y+h/2+1.1,signed(v,2),9.5,anchor='middle',color=WHITE if abs(v)/lim>.72 else INK)
    elif v<0:d.line(x+w/2-.7,y+h/2,x+w/2+.7,y+h/2,INK,.22)

REGIONS={'A':(10,10,156,188),'B':(181,10,128,183),'C':(326,10,184,188),'D':(10,214,146,163),'E':(176,214,220,191),'F':(410,214,100,166),'G1':(10,420,332,145),'G2':(364,420,146,145)}
TITLES={'A':'State identity','B':'PTPN22 recurrence','C':'PTPN22-program architecture','D':'Whole-tumor context','E':'Patient coupling structure','F':'Lineage / technical boundary','G1':'One perturbation experiment','G2':'Human-genetic boundary'}

def build():
    d=D(520,575)
    for tag,(x,y,w,h) in REGIONS.items():d.text(x,y+1,tag,18,True)
    # A: frozen dot area = detection fraction; tone = existing within-marker z.
    lookup={(r['marker'],r['state']):r for r in A};x0,y0,cw,rh=77,47,13.2,4.05
    for j,lab in enumerate(state_labels):d.text(x0+cw*(j+.5)+1,44,lab,10,angle=90)
    block_names={'Progenitor_memory':'Progenitor / memory','Effector_cytotoxic':'Effector / cytotoxic','Checkpoint_exhaustion':'Checkpoint / exhaustion','Proliferation':'Proliferation','Tfh_support':'Tfh support','Regulatory':'Regulatory'}
    blocks={}
    for i,g in enumerate(markers):
        yy=y0+rh*(i+.5);d.text(x0-1,yy+1,g,10,anchor='end')
        block=lookup[g,states[0]]['marker_block'];blocks.setdefault(block,[]).append(i)
        for j,s in enumerate(states):
            r=lookup[g,s];v=num(r,'within_marker_z');f=num(r,'patient_median_positive_fraction');cx=x0+cw*(j+.5)
            assert 0<=f<=1 and abs(v)<=2.5
            radius=1.85*math.sqrt(f);d.dot(cx,yy,radius,color(v,2.5))
            if v<0:d.line(cx+2.4,yy,cx+3.5,yy,GRAY,.2)
            d.mark(panel='A',marker=g,state=s,value=v,detection_fraction=f,n=int(r['n_patients']))
    for name,inds in blocks.items():
        yy=y0+rh*inds[0];end=y0+rh*(inds[-1]+1)
        d.line(57,yy+.6,57,end-.6,GRAY,.3)
        for k,line in enumerate(block_names[name].split(' / ')):d.text(54,(yy+end)/2+(k-.5)*3.6,line,9.5,anchor='end',color=GRAY)
        if inds[0]:d.line(60,yy,156,yy,LIGHT,.2)
    cbar(d,76,181,35,2.5,'Within-marker z')
    for j,f in enumerate([.1,.5,1]):
        cx=124+j*13;d.dot(cx,183,1.85*math.sqrt(f),GRAY);d.text(cx,189,f'{int(f*100)}%',9,anchor='middle')
    d.text(124,194,'Detection',9.5)
    d.text(15,187,'Dash: z < 0',9.5)
    # B: complete six-state strips for both cohorts, no point+CI glyphs.
    blabs=['Effector','Naive / progenitor','Proliferating','Terminal exhaustion','Tfh-like','Treg']
    for j,co in enumerate(['GSE149614','GSE326201']):d.text(248+j*35,29,co,10,True,anchor='middle')
    for i,r in enumerate(B):
        yy=41+i*21;d.text(182,yy+5,blabs[i],10)
        for j,suffix in enumerate(['discovery','replication']):
            xx=233+j*35;v=num(r,'median_patient_rho_'+suffix);lo=num(r,'ci_lower_'+suffix);hi=num(r,'ci_upper_'+suffix)
            cell(d,xx,yy,30,9,v,.15,number=False)
            d.text(xx+15,yy+6,signed(v),11,True,anchor='middle',color=WHITE if abs(v)/.15>.72 else INK)
            d.text(xx+15,yy+13,f'[{lo:+.3f}, {hi:+.3f}]',8.5,anchor='middle')
            d.text(xx+15,yy+18,f"{num(r,'positive_patient_fraction_'+suffix)*100:g}% positive",8.5,anchor='middle',color=GRAY)
            d.mark(panel='B',state=r['state'],cohort=r['dataset_'+suffix],value=v,low=lo,high=hi,positive_fraction=num(r,'positive_patient_fraction_'+suffix),n=int(r['n_patients_'+suffix]))
    cbar(d,234,175,62,.15,'Patient-median Spearman ρ')
    # C: frozen 12 x 6 beta, no clustering and no invented biological families.
    lut={(r['program'],r['state']):r for r in C};xx,yy,cw,rh=395,49,18.5,10.2
    for j,lab in enumerate(state_labels):d.text(xx+cw*(j+.5)+1,46,lab,10,angle=90)
    for i,p in enumerate(programs):
        d.text(xx-3,yy+rh*(i+.5)+1,program_labels[i],10,anchor='end')
        for j,s in enumerate(states):
            r=lut[p,s];v=num(r,'beta_per_1SD_PTPN22');lo=num(r,'ci_lower');hi=num(r,'ci_upper')
            cell(d,xx+j*cw,yy+i*rh,cw,rh,v,.8,outline=r['ci_excludes_zero']=='TRUE')
            d.mark(panel='C',program=p,state=s,value=v,low=lo,high=hi,ci_excludes_zero=r['ci_excludes_zero'],n=int(r['n_patient_states']))
    cbar(d,395,181,48,.8,'β per 1 SD PTPN22')
    d.text(452,184,'Outline: CI excludes 0',9.5)
    if any(not math.isfinite(num(r,'beta_per_1SD_PTPN22')) for r in C):d.text(452,189,'× NA',9.5)
    # D: vertical cohort skyline with slim CI sleeves; frozen meta/PI secondary.
    x0,y0,w,h=29,236,112,81;lo,hi=-1.5,3;fy=lambda v:y0+h*(hi-v)/(hi-lo)
    for v in [-1,0,1,2,3]:
        yy=fy(v);d.text(x0-3,yy+1,str(v),10,anchor='end',color=GRAY)
        d.line(x0,yy,x0+w,yy,GRAY if v==0 else LIGHT,.35 if v==0 else .16)
    d.text(10,228,'Hedges g (R - NR)',10)
    for i,r in enumerate(COHORT):
        cx=x0+14+i*28;v=num(r,'hedges_g');a=num(r,'ci_lower_g');b=num(r,'ci_upper_g')
        d.rect(cx-2,fy(max(0,v)),4,abs(fy(v)-fy(0)),TEAL)
        d.line(cx,fy(a),cx,fy(b),GRAY,.45)
        for end in [a,b]:d.line(cx-2.3,fy(end),cx+2.3,fy(end),GRAY,.32)
        d.text(cx+4,fy(v)-3,signed(v,2),10,True)
        d.text(cx+1,346,r['dataset'],10,angle=90)
        d.mark(panel='D',cohort=r['dataset'],value=v,low=a,high=b,n_R=int(r['n_responder']),n_NR=int(r['n_nonresponder']))
    d.text(10,358,'Pooled',10,True)
    v,a,b=[num(META,k) for k in ['pooled_hedges_g','ci_low','ci_high']];pl,ph=[num(META,k) for k in ['prediction_interval_low','prediction_interval_high']]
    xx,ww=65,80;pos=lambda z:xx+ww*(z+.4)/2.4
    d.line(pos(pl),357,pos(ph),357,LIGHT,2.0);d.line(pos(a),357,pos(b),357,GRAY,.7)
    d.rect(pos(v)-.75,355.5,1.5,3,TEAL);d.line(pos(0),353,pos(0),361,INK,.25)
    d.text(10,367,f'g {v:.3f} [{a:.3f}, {b:.3f}]',10)
    d.text(10,373,f'PI [{pl:.3f}, {ph:.3f}]',10)
    d.mark(panel='D',cohort='PRIMARY_REML_HKSJ',value=v,low=a,high=b,pi_low=pl,pi_high=ph)
    # E: six patient identities in each of three explicit lineage blocks, signed values.
    xx,yy,cw,rh=215,246,28,6.8
    d.text(xx+1.5*cw,230,'Original',12,True,anchor='middle');d.text(xx+4.5*cw+5,230,'Technical-adjusted',12,True,anchor='middle')
    for j,lab in enumerate(['Total','Between','Within']*2):d.text(xx+cw*(j+.5)+(5 if j>=3 else 0),240,lab,10.5,True,anchor='middle')
    elut={(r['lineage'],r['patient'],r['condition'],r['component']):r for r in E}
    for k,lineage in enumerate(['All T cells','CD4','CD8']):
        top=yy+k*(6*rh+5);d.text(176,top+5,lineage,11,True)
        for i,p in enumerate(patient_ids):
            d.text(xx-3,top+rh*(i+.5)+1,p,10,anchor='end')
            for j,component in enumerate(['C_total','C_between','C_within']*2):
                cond='Original' if j<3 else 'Technical-adjusted';r=elut[lineage,p,cond,component];v=num(r,'value')
                cell(d,xx+j*cw+(5 if j>=3 else 0),top+i*rh,cw,rh,v,.5,number=False)
                d.text(xx+cw*(j+.5)+(5 if j>=3 else 0),top+rh*(i+.5)+1,signed(v),10,anchor='middle',color=WHITE if abs(v)/.5>.72 else INK)
                d.mark(panel='E',lineage=lineage,patient=p,condition=cond,component=component,value=v)
    cbar(d,218,389,75,.5,'Signed standardized covariance')
    # F: median total paired bars; exact original summaries, not component sums.
    d.text(410,228,'Median total',11,True)
    for k,(lineage,label) in enumerate([('All T cells','All T: between-dominant'),('CD4','CD4: mixed'),('CD8','CD8: unresolved')]):
        yy=239+k*26;d.text(410,yy,label,10,True)
        for j,cond in enumerate(['PRIMARY','TECHNICAL_ADJUSTED']):
            r=next(r for r in F if r['lineage']==lineage and r['condition']==cond);v=num(r,'median_C_total')
            bar(d,445,yy+6+j*7,47,v,0,.25,h=3.8,col=BLUE if j==0 else TEAL)
            d.text(410,yy+7+j*7,'Original' if j==0 else 'Technical',9.5)
            d.text(507,yy+7+j*7,f'{v:.3f}',9.5,anchor='end')
            d.mark(panel='F',kind='lineage_summary',lineage=lineage,condition=cond,value=v,n=int(r['n']))
    d.text(410,316,'Lineage stability not preserved',9.5)
    d.text(410,327,'State omission: Δmedian between',9.5)
    # All 20 existing single-state removals, same original order, no new deletion analysis.
    for i,r in enumerate(LOSO):
        cx=414+i*4.8;v=num(r,'signed_change_C_between');yy=349-v/.05*13
        d.rect(cx-1, min(349,yy),2,max(abs(349-yy),.1),TEAL if v>=0 else BLUE)
        short=r['removed_state'].split('_')[1].replace('C','')
        d.text(cx,368,short,9,anchor='middle')
        d.mark(panel='F',kind='LOSO',state=r['removed_state'],value=v,final_gate=r['final_gate'],n=int(r['eligible_patient_n']))
    d.line(411,349,510,349,GRAY,.25);d.line(455,332,455,373,LIGHT,.25)
    d.text(433,378,'CD4 removals',9.5,anchor='middle');d.text(483,378,'CD8 removals',9.5,anchor='middle')
    d.text(410,337,'+0.05',8.5,color=GRAY);d.text(410,363,'-0.05',8.5,color=GRAY)
    # G1: all 35 genes x 3 conditions from one experiment, frozen orientation.
    gstart=len(d.ops)
    d.text(28,401,'GSE314342',11,True)
    genes=[r['gene'] for r in sorted([r for r in G1 if r['condition']=='Rest'],key=lambda r:int(r['gene_rank']))]
    conditions=['Rest','Stim8hr','Stim48hr'];glut={(r['gene'],r['condition']):r for r in G1}
    xx,yy,cw,rh=43,429,8.25,10
    for j,g in enumerate(genes):d.text(xx+cw*(j+.5)+1,426,g,10,angle=90)
    lim=2.5
    assert all(abs(num(r,'oriented_effect'))<=lim for r in G1)
    for i,c in enumerate(conditions):
        d.text(10,yy+rh*(i+.5)+1,{'Rest':'Rest','Stim8hr':'8 h','Stim48hr':'48 h'}[c],10.5,True)
        for j,g in enumerate(genes):
            r=glut[g,c];v=num(r,'oriented_effect');cell(d,xx+j*cw,yy+i*rh,cw,rh,v,lim,number=False)
            d.mark(panel='G1',kind='gene_effect',gene=g,condition=c,value=v,native_effect=num(r,'native_effect'),orientation_factor=num(r,'orientation_factor'))
    cbar(d,43,465,62,lim,'Oriented effect: -log_fc (CRISPRi - NTC)')
    d.text(10,489,'Program shift',10.5,True);d.text(235,489,'Mean oriented log_fc',10.5)
    xx,ww,lo,hi=61,143,-.25,.25
    for i,r in enumerate(GSH):
        ref=next(g for g in GR if g['reference']==r['system']);yy=500+i*10;v=num(r,'AF35_PROGRAM_SHIFT');a=num(ref,'random_q025');b=num(ref,'random_q975')
        pos=lambda z:xx+ww*(z-lo)/(hi-lo)
        d.rect(pos(a),yy-2,pos(b)-pos(a),4,'#E5EAEC');bar(d,xx,yy,ww,v,lo,hi,h=2.5,col=TEAL if v>=0 else BLUE)
        d.line(pos(0),yy-4,pos(0),yy+4,GRAY,.25)
        d.text(10,yy+1,{'Rest':'Rest','Stim8hr':'8 h','Stim48hr':'48 h'}[r['condition']],10)
        d.text(210,yy+1,signed(v),10);d.text(274,yy+1,signed(num(r,'mean_raw_oriented_effect')),10,anchor='middle')
        d.mark(panel='G1',kind='program_shift',condition=r['condition'],value=v,reference_low=a,reference_high=b,raw_mean=num(r,'mean_raw_oriented_effect'))
    haxis(d,xx,527,ww,lo,hi,[-.2,0,.2]);d.text(10,544,'Shading: random-gene 2.5-97.5% range',10)
    # G2: every donor, separate group strips; no genotype recoding or carrier language.
    d.text(386,401,'GSE94859',11,True);d.text(364,413,'AF35 Δ (Activated - Naive; z units)',10.5)
    xx,y0,ww,hh,lo,hi=389,425,109,74,-3.5,1.0;fy=lambda v:y0+hh*(hi-v)/(hi-lo)
    for v in [-3,-2,-1,0,1]:
        yy=fy(v);d.text(xx-4,yy+1,str(v),10,anchor='end');d.line(xx,yy,xx+ww,yy,GRAY if v==0 else LIGHT,.3 if v==0 else .16)
    for j,gen in enumerate(['CC','TT']):
        selected=[r for r in G2 if r['genotype']==gen]
        for i,r in enumerate(selected):
            cx=xx+10+j*57+(i%4)*7;cy=fy(num(r,'delta'))
            # Deterministic donor lanes use no random jitter and do not alter y values.
            d.dot(cx,cy,1.4,INK,open=gen=='TT')
            d.text(cx+2,cy+1,r['individual_id'].replace('individual',''),8.5,color=GRAY)
            d.mark(panel='G2',kind='donor_delta',donor=r['individual_id'],genotype=gen,value=num(r,'delta'))
        d.text(xx+22+j*57,506,gen,11,True,anchor='middle')
    d.text(364,516,'TT - CC',10.5,True);d.text(454,516,'95% donor CI',10)
    for j,(label,vk,lk,hk,lo,hi) in enumerate([('D','D_TT_minus_CC','bootstrap_D_CI95_low','bootstrap_D_CI95_high',-2,1),('g','Hedges_g','bootstrap_g_CI95_low','bootstrap_g_CI95_high',-2,1)]):
        yy=526+j*12;v,a,b=[num(GS,k) for k in [vk,lk,hk]];x,w=387,63;pos=lambda z:x+w*(z-lo)/(hi-lo)
        d.text(364,yy+1,label,10.5);d.line(pos(0),yy-3,pos(0),yy+3,LIGHT,.25)
        d.line(pos(a),yy,pos(b),yy,GRAY,.35);d.rect(pos(v)-.7,yy-1.2,1.4,2.4,GRAY)
        for z in [a,b]:d.line(pos(z),yy-.8,pos(z),yy+.8,GRAY,.25)
        d.text(510,yy+1,f'{v:.3f} [{a:.3f}, {b:.3f}]',9,anchor='end')
        d.mark(panel='G2',kind='genotype_interval',field=vk,value=v,low=a,high=b)
    for i in range(gstart,len(d.ops)):
        op=list(d.ops[i]);op[2]+=20
        if op[0]=='line':op[4]+=20
        d.ops[i]=tuple(op)
    return d

if __name__=='__main__':
    d=build();d.render(OUT/(PREFIX+'_ASSEMBLY_READY'))
    # Same geometry/scientific objects; author-review labels occupy only reserved outer header space.
    review=D(d.w,d.h);review.ops=list(d.ops);review.marks=list(d.marks)
    for tag,(x,y,w,h) in REGIONS.items():
        if tag not in ['G1','G2']:review.text(x+14,y+1,TITLES[tag],11,True,color=GRAY)
    review.render(OUT/(PREFIX+'_AUTHOR_REVIEW'))
    dump(OUT/'FIGURE_DIMENSIONS_AND_FONTS.json',dict(width_mm=d.w,height_mm=d.h,regions=REGIONS,fonts={k:str(v) for k,v in base.FONTFILES.items()},minimum_font_pt=min(o[4] for o in d.ops if o[0]=='text'),panel_tag_font_pt=18,scientific_marks=len(d.marks),new_analyses=0))
    doc=fitz.open(OUT/(PREFIX+'_ASSEMBLY_READY.pdf'));doc[0].get_pixmap(dpi=100,alpha=False).save(QA/'NATIVE_OVERVIEW.png')
    print(f'RENDERED_ONE_MASTER; marks={len(d.marks)}; dimensions={d.w}x{d.h} mm')
