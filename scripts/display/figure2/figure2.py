from common import *
from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap

q=Work(2,560,720)
genes=[r['gene'] for r in q.src('AF35_ORDERED_GENES.tsv')]
ann=q.src('AF35_FUNCTIONAL_ROLE_ANNOTATION.tsv');coupling=q.src('AF35_PTPN22_GENE_COUPLING.tsv')
assert len(genes)==35 and [r['gene'] for r in ann]==genes and 'PTPN22' not in genes
cohorts=['ERP117672','GSE302495']
mat={co:np.array([[float(r[g]) for g in genes] for r in q.src('AF35_PROGRAM_COHERENCE_GENE_CORRELATION_'+co+'.tsv')]) for co in cohorts}
pairs=q.src('AF35_PROGRAM_COHERENCE_COVARIANCE_VECTORS.tsv')
cross=q.src('AF35_PROGRAM_COHERENCE_CROSS_COHORT.tsv')[0]
hs=q.src('AF35_PROGRAM_COHERENCE_WITHIN_COHORT_SUMMARY.tsv')
splits=q.src('AF35_PROGRAM_COHERENCE_SPLIT_HALF.tsv')
single=q.src('F2E_SINGLE_DELETION_SCORE_ONLY.tsv')
multi=q.src('AF35_GENE_NUMBER_SUMMARY.tsv');draws=q.src('AF35_GENE_NUMBER_SCORE_STABILITY.tsv')
random=q.src('AF35_PROGRAM_COHERENCE_RANDOM_REFERENCE_RESULTS.tsv');ref=q.src('AF35_PROGRAM_COHERENCE_RANDOM_REFERENCE_SUMMARY.tsv')
assert len(pairs)==595 and len(single)==70 and len(random)==2000
assert [(r['gene_i'],r['gene_j']) for r in pairs]==[(genes[i],genes[j]) for i in range(35) for j in range(i+1,35)]
for r in pairs:
    for co in cohorts:assert float(r[co])==mat[co][genes.index(r['gene_i']),genes.index(r['gene_j'])]

# A: source annotation layers are membership flags, not new biological grouping.
q.panel('A',12,13,536,92,'AF35 members and PTPN22 coupling')
layers=list(dict.fromkeys(z for r in ann for z in r['interpretive_layers'].split(';')))
ax=q.native_ax(31,36,40,181)
flags=np.array([[int(l in r['interpretive_layers'].split(';')) for r in ann] for l in layers])
layer_colors=['#C84F55','#D68A32','#2B90A3','#715BAD','#2A6AAF','#308D6B','#B34587','#606C84']
for k,l in enumerate(layers):
    for j,g in enumerate(genes):
        if flags[k,j]:ax.add_patch(Rectangle((k-.34,j-.32),.68,.64,facecolor=layer_colors[k],edgecolor='none'))
ax.set(xlim=(-.5,len(layers)-.5),ylim=(34.5,-.5));ax.set_yticks(range(35),genes)
ax.set_xticks(range(len(layers)),[{'Negative':'Feedback','NFkB':'NFκB'}.get(l,l) for l in layers],rotation=90);ax.xaxis.tick_top();ax.tick_params(length=0)
ax.spines[['left','bottom']].set_visible(False)
for i,l in enumerate(layers):
    for j,g in enumerate(genes):q.mark('A',kind='annotation',gene=g,layer=l,present=int(flags[i,j]))
ax=q.native_ax(78,36,27,181)
vals=np.array([[float(next(r for r in coupling if r['gene']==g and r['cohort']==co)['PTPN22_spearman']) for g in genes] for co in cohorts])
im=ax.pcolormesh(np.arange(3)-.5,np.arange(36)-.5,vals.T,cmap=DIV,vmin=-1,vmax=1);ax.invert_yaxis();ax.set_xticks([0,1],['ERP','GSE']);ax.xaxis.tick_top();ax.set_yticks([]);ax.tick_params(length=0);ax.spines[['left','bottom']].set_visible(False)
for i in range(2):
    for j in range(35):
        ax.text(i,j,f'{vals[i,j]:+.2f}',ha='center',va='center',fontsize=8,color='white' if abs(vals[i,j])>.7 else INK)
        q.mark('A',kind='coupling',gene=genes[j],cohort=cohorts[i],value=vals[i,j])
q.native_text(78,226,'PTPN22 ρ',9)
q.native_text(31,226,'35 members',9)

# B: both full symmetric matrices represented once each in their triangle.
q.panel('B',12,124,286,289,'AF35 internal structure')
ax=q.ax(46,151,244,244)
dual=np.where(np.triu(np.ones((35,35)),1).astype(bool),mat[cohorts[0]],mat[cohorts[1]])
im=ax.pcolormesh(np.arange(36)-.5,np.arange(36)-.5,dual,cmap=DIV,vmin=-1,vmax=1);ax.invert_yaxis();ax.set_xticks(range(35),genes,rotation=90);ax.xaxis.tick_top();ax.set_yticks(range(35),genes);ax.tick_params(length=0,pad=3);ax.spines[['left','bottom']].set_visible(False)
for i in range(35):
    for j in range(35):
        co=cohorts[0] if j>i else cohorts[1]
        if i!=j and dual[i,j]<0:ax.plot([j-.16,j+.16],[i,i],color=INK,lw=.7)
        q.mark('B',gene_i=genes[i],gene_j=genes[j],cohort=co,value=float(dual[i,j]))
q.text(48,404,'Upper: ERP117672    Lower: GSE302495',10)
cb=q.fig.colorbar(im,cax=q.ax(47,408,80,2),orientation='horizontal',ticks=[-1,0,1]);cb.set_label('Spearman ρ',labelpad=2)
q.text(164,413,'Dash: ρ < 0',9)

q.panel('C',318,124,230,193,'Cross-cohort pair conservation')
ax=q.ax(344,145,180,155)
x=np.array([float(r[cohorts[0]]) for r in pairs]);y=np.array([float(r[cohorts[1]]) for r in pairs])
hexes=ax.hexbin(x,y,gridsize=26,extent=(-1,1,-1,1),mincnt=1,cmap=SEQ,linewidths=0)
low=np.floor(min(x.min(),y.min())*10)/10
ax.plot([-1,1],[-1,1],ls='--',lw=.7,color=GRAY);ax.set(xlim=(low,1),ylim=(low,1),xlabel='ERP117672 gene-pair ρ',ylabel='GSE302495 gene-pair ρ')
cb=q.fig.colorbar(hexes,cax=q.ax(529,145,3,155));cb.set_label('Gene-pair count')
q.text(345,137,f"r = {float(cross['covariance_pearson']):.3f}    ρ = {float(cross['covariance_spearman']):.3f}    595 pairs",10)
assert int(hexes.get_array().sum())==595
for r in pairs:q.mark('C',**r)

q.panel('D',318,330,230,87,'Measurement reliability')
metrics=[('median_pairwise','Median pairwise ρ'),('PC1_variance','PC1 variance'),('split_half_median','Split-half ρ'),('spearman_brown_median','Spearman–Brown')]
for i,(key,label) in enumerate(metrics):
    ax=q.ax(405,345+i*16,120,11)
    for j,r in enumerate(hs):
        v=float(r[key]);ax.barh(1-j,v,height=.32,color=[BLUE,TEAL][j]);ax.text(v+.015,1-j,f'{v:.3f}',va='center',fontsize=9)
        q.mark('D',cohort=r['cohort'],metric=key,value=v)
    ax.set(xlim=(0,1.13),ylim=(-.5,1.5));ax.set_yticks([1,0],['ERP','GSE']);ax.set_xticks([]);ax.spines[['left','bottom']].set_visible(False);ax.tick_params(length=0)
    q.text(320,351+i*16,label,10)

q.panel('E',12,441,288,117,'Single-member score ablation')
for i,co in enumerate(cohorts):
    rs=sorted([r for r in single if r['cohort']==co],key=lambda r:int(r['gene_rank']));assert [r['omitted_gene'] for r in rs]==genes
    vs=np.array([float(r['score_rho']) for r in rs]);ax=q.ax(44,453+i*39,248,30)
    assert vs.min()>=.94
    ax.bar(np.arange(35),vs,width=.67,color=[BLUE,TEAL][i]);ax.set(ylim=(.94,1.002),xlim=(-.6,34.6),ylabel='Score-to-full ρ')
    ax.set_yticks([.94,.97,1]);ax.set_xticks(range(35),genes if i else ['']*35,rotation=90);ax.tick_params(axis='x',length=0);ax.text(.01,.92,co,transform=ax.transAxes,fontsize=10)
    for r in rs:q.mark('E',**r)

q.panel('F',318,441,230,117,'Retained-member robustness')
for i,co in enumerate(cohorts):
    rs=sorted([r for r in multi if r['cohort']==co and int(r['k'])<35],key=lambda r:int(r['k']))
    assert [int(r['k']) for r in rs]==[5,10,15,20,25,30]
    ax=q.ax(346+i*98,465,83,70)
    for r in rs:
        k=int(r['k']);v,a,b=map(float,[r['rho_median'],r['rho_q025'],r['rho_q975']])
        assert len([z for z in draws if z['cohort']==co and int(z['k'])==k])==2000
        ax.add_patch(Rectangle((k-1.8,a),3.6,b-a,facecolor=[BLUE,TEAL][i],alpha=.22,edgecolor='none'))
        ax.plot([k-1.8,k+1.8],[v,v],color=[BLUE,TEAL][i],lw=2)
        q.mark('F',cohort=co,k=k,median=v,q025=a,q975=b,draws=2000)
    floor=min(float(r['rho_q025']) for r in multi if int(r['k'])<35)
    full=next(r for r in multi if r['cohort']==co and int(r['k'])==35)
    ax.plot([33.5,36.5],[float(full['rho_median'])]*2,color=INK,lw=2)
    q.mark('F',cohort=co,k=35,median=float(full['rho_median']),draws=1)
    ax.set(xlim=(2,38),ylim=(np.floor(floor*10)/10,1.01),xlabel='Retained genes',ylabel='Score-to-full ρ' if i==0 else '');ax.set_xticks([5,10,15,20,25,30,35]);ax.set_title(co,fontsize=10)

q.panel('G',12,580,536,129,'Matched-random structural reference')
short={'median_pairwise':'Pairwise ρ','PC1_variance':'PC1 variance','median_gene_to_rest':'Gene-to-rest ρ','covariance_pearson':'Pair conservation r','covariance_spearman':'Pair conservation ρ','connectivity_pearson':'Connectivity r','connectivity_spearman':'Connectivity ρ'}
for i,r in enumerate(ref):
    key=r['metric'];ax=q.ax(35+(i%5)*104,596+(i//5)*59,84,38)
    v=np.array([float(z[key]) for z in random]);obs=float(r['AF35_observed'])
    ax.hist(v,bins=30,color='#AEBEC0',edgecolor='white',linewidth=.2);ax.axvline(obs,color=TEAL,lw=1.5)
    ax.text(obs,ax.get_ylim()[1]*.9,f'{obs:.3f}',ha='right' if obs>.9 else 'left',fontsize=9,color=TEAL)
    label=('ERP ' if key.startswith('ERP') else 'GSE ' if key.startswith('GSE') else '')+short[key.split('__')[-1]]
    ax.set_xlabel(label);ax.set_ylabel('Programs' if i%5==0 else '');ax.tick_params(axis='y',labelsize=9)
    for z in random:q.mark('G',metric=key,program=z['program_id'],value=float(z[key]))
    q.mark('G',metric=key,observed=obs,n_ge=int(r['n_ge_AF35']))
q.save()
