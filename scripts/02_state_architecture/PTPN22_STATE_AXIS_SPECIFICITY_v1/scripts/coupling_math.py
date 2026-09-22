import numpy as np
from scipy.stats import rankdata
def correlation(x,y,spearman=False):
    x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float)
    if spearman:x=rankdata(x,method='average');y=rankdata(y,method='average')
    x=x-x.mean();y=y-y.mean();den=np.sqrt(np.sum(x*x)*np.sum(y*y))
    return np.sum(x*y)/den if den>0 else np.nan
def decompose(x,y,states):
    x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float);states=np.asarray(states)
    if len(x)<2 or not np.isfinite(x).all() or not np.isfinite(y).all() or np.var(x)==0 or np.var(y)==0:return None
    xz=(x-x.mean())/x.std(ddof=0);yz=(y-y.mean())/y.std(ddof=0)
    total=np.mean((xz-xz.mean())*(yz-yz.mean()));between=within=0.;rx=xz.copy();ry=yz.copy();detail=[]
    for state in np.unique(states):
        m=states==state;n=int(m.sum());weight=n/len(x);mx=xz[m].mean();my=yz[m].mean()
        cv=np.mean((xz[m]-mx)*(yz[m]-my))
        cb=weight*(mx-xz.mean())*(my-yz.mean());cw=weight*cv
        between+=cb;within+=cw;rx[m]-=mx;ry[m]-=my
        detail.append(dict(state=str(state),cell_n=n,weight=weight,mean_PTPN22=float(x[m].mean()),mean_AF35=float(y[m].mean()),mean_Xz=mx,mean_Yz=my,Cov_state=cv,between_contribution=cb,within_contribution=cw))
    error=total-between-within
    assert abs(error)<1e-12
    return dict(cell_n=len(x),state_n=len(detail),PTPN22_variance=float(x.var()),AF35_variance=float(y.var()),C_total=total,C_between=between,C_within=within,identity_error=error,residual_Pearson=correlation(rx,ry),residual_Spearman=correlation(rx,ry,True),states=detail)
def technical_residual(x,y,umi,genes):
    design=np.column_stack([np.ones(len(x)),np.log1p(umi),np.log1p(genes)])
    centered=design[:,1:]-design[:,1:].mean(axis=0)
    sd=centered.std(axis=0);nonzero=sd>0
    design=np.column_stack([np.ones(len(x)),centered[:,nonzero]/sd[nonzero]])
    values=np.column_stack([x,y]);coef=np.linalg.lstsq(design,values,rcond=None)[0]
    res=values-design@coef
    return res[:,0],res[:,1],int(np.linalg.matrix_rank(design))
