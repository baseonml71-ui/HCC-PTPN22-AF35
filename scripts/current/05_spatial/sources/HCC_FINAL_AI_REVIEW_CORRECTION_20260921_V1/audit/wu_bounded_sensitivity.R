invisible(Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8'))
options(stringsAsFactors=FALSE,digits=16)
suppressPackageStartupMessages({library(data.table);library(Matrix);library(FNN);library(rhdf5)})
O <- 'HCC_FINAL_AI_REVIEW_CORRECTION_20260921_V1'
src <- 'HCC_ICI_project/04_scripts/spatial_specificity_rescue/01_spatial_specificity_rescue.R'
for(e in parse(src)) if(is.call(e)&&identical(e[[1]],as.name('<-'))&&as.character(e[[2]]) %in% c('z_safe','conditional_residuals','moran_matrix','csr_library_sizes','group_effect')) eval(e)
p <- 'HCC_ICI_project/03_processed_data/spatial_phase8/Visium-ST/visium_all.h5ad'
dec<-function(k){a<-h5read(p,paste0('/obs/',k,'/categories'));b<-h5read(p,paste0('/obs/',k,'/codes'));a[b+1L]}
obs<-data.table(sample_id=as.character(dec('sample_id')),state=as.character(dec('diagnosis')),response=as.character(dec('Response')),in_tissue=as.character(dec('in_tissue')))
obs[,response:=ifelse(response=='Responder','Responder','Nonresponder')]
coords<-t(h5read(p,'/obsm/spatial'));storage.mode(coords)<-'double'
ptr<-h5read(p,'/X/indptr'); genes<-as.character(h5read(p,'/var/_index'))
af<-fread('AF35_RESPONSE_BLIND_PROGRAM_COHERENCE_v1/inputs/AF35_ORDERED_GENES.tsv')$gene
co<-fread('HCC_ICI_project/01_metadata/PHASE3_PTPN22_SPECIFICITY_COMPARATOR_PROGRAMS.tsv')
broad<-intersect(co[program=='broad_immune_abundance' & !startsWith(gene,'DERIVED:'),gene],genes)
stopifnot(length(af)==35L,!('PTPN22'%in%af))
frozen<-fread('HCC_ICI_project/05_results/tables/spatial_specificity_rescue/PTPN22_CONDITIONAL_SPATIAL_EFFECTS.tsv')[cohort=='Wu_Zenodo_19123188' & record_type=='PATIENT']
aud<-list();eff<-list()
for(s in unique(obs[state=='Pre',sample_id])){
 ii<-which(obs$sample_id==s);stopifnot(all(diff(ii)==1));a<-ptr[min(ii)]+1;b<-ptr[max(ii)+1]
 ix<-as.integer(h5read(p,'/X/indices',index=list(a:b)));v<-as.numeric(h5read(p,'/X/data',index=list(a:b)))
 row<-rep.int(seq_along(ii),diff(ptr[min(ii):(max(ii)+1)])); tot<-csr_library_sizes(v,row,length(ii));det<-tabulate(row,nbins=length(ii))
 selected<-(ix+1L) %in% match(c(af,broad),genes); lv<-log1p(v*10000/tot[row])
 mat<-sparseMatrix(i=row[selected],j=ix[selected]+1L,x=lv[selected],dims=c(length(ii),length(genes)))
 target<-as.numeric(rowMeans(mat[,match(af,genes)]));imm<-as.numeric(rowMeans(mat[,match(broad,genes)]))
 keep<-obs$in_tissue[ii]=='1' & tot>0
 getI<-function(k){r<-conditional_residuals(matrix(target[k],ncol=1),imm[k],log1p(tot[k]),log1p(det[k]))$residual;nn<-FNN::get.knn(coords[ii[k],,drop=FALSE],k=6L)$nn.index;as.numeric(moran_matrix(r,nn))}
 orig<-getI(rep(TRUE,length(ii)));fv<-frozen[sample_id==s,conditional_moran];stopifnot(length(fv)==1,abs(orig-fv)<1e-10);sen<-getI(keep)
 aud[[s]]<-data.table(sample_id=s,patient_id=sub('cytassist_([0-9]+)_pre','\\1',s),response=obs$response[ii[1]],total_source_spots=length(ii),tissue_included_spots=sum(obs$in_tissue[ii]=='1'),zero_total_count_spots=sum(tot==0),zero_count_spots_in_original_analysis=sum(tot==0),spots_entering_original_residualization=length(ii),spots_entering_original_knn=length(ii),out_of_tissue_spots=sum(obs$in_tissue[ii]!='1'),explicit_analysis_tissue_filter='NONE; source object in_tissue=1 for every retained row',coordinate_source='visium_all.h5ad /obsm/spatial',disconnected_tissue_islands='NOT_RECORDED',neighbor_distance_summary='NOT_RECORDED',classification=if(any(!keep))'C. MATERIAL_INPUT_RISK_REQUIRES_BOUNDED_SENSITIVITY' else 'B. SOURCE_FILTER_CONFIRMED',sensitivity_retained_spots=sum(keep))
 eff[[s]]<-data.table(sample_id=s,patient_id=aud[[s]]$patient_id,response=obs$response[ii[1]],frozen_moran=fv,recovered_original_moran=orig,filtered_moran=sen,change=sen-fv)
 cat(s,'zero=',sum(tot==0),'original=',orig,'filtered=',sen,'\n')
}
audit<-rbindlist(aud); effects<-rbindlist(eff);stopifnot(sum(audit$zero_total_count_spots)==397)
fwrite(audit,file.path(O,'WU_PRETREATMENT_SPATIAL_INPUT_AUDIT.tsv'),sep='\t')
fwrite(effects,file.path(O,'recovered_sources/WU_ZERO_COUNT_SENSITIVITY_PATIENT_METRICS.tsv'),sep='\t')
n_boot<-20000L;set.seed(20260921L);s<-group_effect(effects$filtered_moran,effects$response);s[,analysis:='single_source_supported_zero_count_exclusion'];s[,seed:=20260921L]
fwrite(s,file.path(O,'recovered_sources/WU_ZERO_COUNT_SENSITIVITY_COHORT_EFFECT.tsv'),sep='\t');print(s)
capture.output(sessionInfo(),file=file.path(O,'audit/WU_SENSITIVITY_SESSION_INFO.txt'))
