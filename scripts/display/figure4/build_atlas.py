from pathlib import Path
import sys,json,hashlib
sys.path.insert(0,'绘图视觉/PARALLEL_FINAL_FIGURE_PRODUCTION_V1/_runtime')
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.patches import Rectangle,Ellipse
from matplotlib.font_manager import FontProperties
O=Path(__file__).parent;R=O.parent.parent;cfg=json.loads((O/'geometry_style.json').read_text());used=[]
def read(p):
 used.append({'name':p.name,'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()});return pd.read_csv(p,sep='\t')
base=R/'HCC_ICI_MAIN_FIGURE_FROZEN_SOURCE_PAYLOAD_V2/Figure5'
table=read(base/'display_support/F5C_ALL_PATIENTS_AND_FROZEN_EFFECTS.tsv');patients=table[table.record_type=='PATIENT'];cohorts=list(patients.cohort.unique())
spots=[read(R/'绘图视觉/PARALLEL_FINAL_FIGURE_PRODUCTION_V1/display_support/F5C_GSE_FROZEN_RESIDUAL_SPOTS.tsv'),read(base/'WU_SAVED_RESIDUAL_SPOTS_EXPORT.tsv')]
wp=R/'.codex_tmp/figure5_pilot2/display_parameters.json';windows={x['sample']:x for x in json.loads(wp.read_text()) if 'sample' in x};used.append({'name':wp.name,'path':str(wp),'sha256':hashlib.sha256(wp.read_bytes()).hexdigest()})
gpath=R/'HCC_ICI_project/05_results/tables/phase2/PTPN22_SPATIAL_PATIENT_METRICS.tsv';wpath=R/'HCC_ICI_project/05_results/tables/spatial_phase8/WU_VISIUM_PATIENT_TIMEPOINT_METRICS.tsv'
assert hashlib.sha256(gpath.read_bytes()).hexdigest()=='1d5088d976fc9b1533d72fdeb2d49949581c4d38098a2c4cea74e0f641c856c9'
assert hashlib.sha256(wpath.read_bytes()).hexdigest()=='30ba771b194706022385db2207b07eec1f61c314b2a7006e2e345c7196ee7af4'
gs=read(gpath).query("metric == 'state_abundance_mean_log1p_CPM' and niche == 'PTPN22_associated_activation_feedback'").set_index('sample_ID');ws=read(wpath).query("state == 'Pre'").set_index('sample_id');levels={}
for _,r in patients.iterrows():
 src=gs.loc[r.sample_id] if r.cohort==cohorts[0] else ws.loc[r.sample_id];assert int(src.n_spots)==int(r.n_spots);assert (src.response=='Responder')==(r.response=='Responder');levels[r.sample_id]=float(src['value'] if r.cohort==cohorts[0] else src.primary_mean_log1p_cp10k)
assert len(levels)==18
plt.rcParams.update({'font.family':'Arial','font.size':14,'text.color':'black','pdf.fonttype':42,'svg.fonttype':'none'})
W,H=cfg['atlas_canvas_mm'];fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white');audit=[]
signed=LinearSegmentedColormap.from_list('conditional',['#1555D4','#F2EEE9','#EF493D']);lm=LinearSegmentedColormap.from_list('level',['#F3FAFE','#59B6DA']);im=LinearSegmentedColormap.from_list('moran',['#F6F2FC','#AD89DC']);ln=Normalize(0,max(levels.values()));inn=Normalize(0,float(patients.conditional_moran.max()))
def ax(x,y,w,h):return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def text(x,y,s,size=15,bold=False,ha='left'):return fig.text(x/W,1-y/H,s,fontsize=size,fontweight='bold' if bold else 'normal',ha=ha,va='center',color='black')
def patient_label(x,y,t,response):
 font=FontProperties(fname='fonts/seguisb.ttf',size=17);fig.text(x/W,1-y/H,t,fontproperties=font,ha='center',va='center')
 width=fig.canvas.get_renderer().get_text_width_height_descent(t,font,False)[0]/fig.dpi*25.4
 fig.add_artist(Ellipse(((x-width/2-2.4)/W,1-y/H),1.7/W,1.7/H,transform=fig.transFigure,fc='#222222' if response=='Responder' else 'white',ec='#222222',lw=.6))
def colorbar(x,y,lim):
 cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(-lim,lim),cmap=signed),cax=ax(x,y,64,2.8),orientation='horizontal',extend='both');cb.solids.set_rasterized(False);cb.solids.set_edgecolor('face');cb.solids.set_linewidth(.25);cb.set_ticks([-lim,0,lim]);cb.set_ticklabels([f'−{lim:.2f}','0',f'+{lim:.2f}']);cb.ax.tick_params(labelsize=12,length=2);cb.outline.set_linewidth(.6);text(x+32,y-5,'Residual SD',13,ha='center')
text(8,8,'D',28,True)
for ci,co in enumerate(cohorts):
 hy=8 if ci==0 else 143;text(26,hy,'GSE238264 · post-treatment · 7 patients' if ci==0 else 'Wu · pretreatment · 11 patients',20,True)
 r=table[table.record_type=='COHORT_SUMMARY'].iloc[ci];text(270,hy+12,f'R − NR = {r.conditional_R_minus_NR:+.3f}   95% CI [{r.bootstrap_ci_low:+.3f}, {r.bootstrap_ci_high:+.3f}]   P = {r.exact_permutation_p_two_sided:.4f}',15)
 colorbar(524,hy+4,windows['HCC1R']['limits'][1] if ci==0 else 4.5)
 for i,(_,r) in enumerate(patients[patients.cohort==co].iterrows()):
  p=windows[r.sample_id];d=spots[ci][spots[ci].sample_id==r.sample_id];assert len(d)==p['spots']==int(r.n_spots)
  x0,y0,x1,y1=p['window']
  if ci==0:
   cellw=584/7;xx=8+i*cellw;yy=43;bw=cellw-5;bh=80;ratio=min(bw/(x1-x0),bh/(y1-y0))/min((384/7-2)/(x1-x0),50/(y1-y0))
  else:
   row,col=divmod(i,6);cellw=584/6;bw=(384/11-2)*1.25*1.12*1.6;bh=43*1.25*1.12*1.6;xx=8+col*cellw+(cellw-bw)/2;yy=178 if row==0 else 310;ratio=1.4*1.6
  fig.add_artist(Rectangle(((xx-2)/W,1-(yy+bh+2)/H),(bw+4)/W,(bh+13)/H,transform=fig.transFigure,fc='#F7F8F9',ec='none',zorder=-1))
  a=ax(xx,yy,bw,bh);a.set(xlim=(x0,x1),ylim=(y1,y0));a.set_aspect('equal');a.axis('off');a.scatter(d.x,d.y,c=d.conditional_score,cmap=signed,norm=Normalize(*p['limits']),s=(p['diameter_pt']*ratio)**2,alpha=1,linewidths=0,edgecolors='none')
  lab=r.sample_id+' · '+('R' if r.response=='Responder' else 'NR') if ci==0 else r.sample_id.split('_')[1]+' · '+('R' if r.response=='Responder' else 'NR');patient_label(xx+bw/2,yy-6,lab,r.response)
  actual_scale=min(bw/(x1-x0),bh/(y1-y0));audit.append({'sample':r.sample_id,'response':r.response,'cohort':co,'Level':levels[r.sample_id],'Moran_I':float(r.conditional_moran),'points':len(d),'NA':int(d.conditional_score.isna().sum()),'source_window':p['window'],'limits':p['limits'],'display_box_mm':[xx,yy,bw,bh],'scale_x_mm_per_source':actual_scale,'scale_y_mm_per_source':actual_scale,'spot_size_scale_vs_v1':ratio,'rotation':0,'y_down':True})
# One shared summary per cohort, exact source order matching the left atlas.
for ci,co in enumerate(cohorts):
 sy=28 if ci==0 else 163
 text(615,sy,'Patient',15,True);text(678,sy,'Level',15,True,ha='center');text(731,sy,'Moran I',15,True,ha='center')
 for j,(_,r) in enumerate(patients[patients.cohort==co].iterrows()):
  yy=sy+13+j*(12 if ci==0 else 18)
  label=r.sample_id if ci==0 else r.sample_id.split('_')[1]+' · '+('R' if r.response=='Responder' else 'NR')
  text(615,yy,label,14)
  for x,val,cmap,norm in [(655,levels[r.sample_id],lm,ln),(708,float(r.conditional_moran),im,inn)]:
   fig.add_artist(Rectangle((x/W,1-(yy+5)/H),46/W,10/H,transform=fig.transFigure,fc=cmap(norm(val)),ec='none'))
   text(x+23,yy,f'{val:.3f}',14,ha='center')
for y in [131]:fig.add_artist(plt.Line2D([8/W,754/W],[1-y/H]*2,transform=fig.transFigure,color='#DCE1E4',lw=.6))
for x,cmap,norm,title in [(620,lm,ln,'Level'),(686,im,inn,'Moran I')]:
 cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=ax(x,414,52,2.2),orientation='horizontal');cb.solids.set_rasterized(False);cb.solids.set_edgecolor('face');cb.solids.set_linewidth(.25);cb.set_ticks([0,norm.vmax]);cb.set_ticklabels(['0',f'{norm.vmax:.2f}']);cb.ax.tick_params(labelsize=11,length=2);cb.outline.set_visible(False);text(x+26,410,title,12,ha='center')
for ext in ['svg','pdf','png']:fig.savefig(O/('FIGURE4_PATIENT_ATLAS_PLAN2_V1.'+ext),dpi=cfg['core_preview_dpi'])
plt.close(fig)
(O/'ATLAS_SOURCE_BINDINGS.json').write_text(json.dumps(used,indent=2),encoding='utf-8');(O/'ATLAS_GEOMETRY.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
pd.DataFrame(audit).to_csv(O/'PATIENT_ATTACHED_READOUTS.tsv',sep='\t',index=False)
print('Atlas: 7 + 11 frozen patients, attached source Level/Moran, unchanged windows and color limits.')
