options(stringsAsFactors=FALSE,digits=16)
suppressPackageStartupMessages({library(data.table);library(digest);library(jsonlite);library(logistf)})
out <- 'CLINICAL_ANCHOR_READOUT_DECOMPOSITION_v1'; p <- 'HCC_ICI_project'
put <- function(x,n) fwrite(x,file.path(out,n),sep='\t',na='NA')
lock<-fromJSON(file.path(out,'MOLECULAR_GEOMETRY_LOCK.json'))
for(n in names(lock$files)) stopifnot(digest(file.path(out,n),algo='sha256',file=TRUE)==lock$files[[n]])
# Evaluate only frozen pure function definitions; never source the old mutating workflow.
source_script <- file.path(p,'04_scripts/phase15c_public_benchmark/01_run_phase15c_benchmark.R')
for(e in parse(source_script)) if(is.call(e)&&identical(e[[1]],as.name('<-'))&&as.character(e[[2]])%in%c('hedges_g','safe_z','fit_univariate_firth')) eval(e)
erpm <- fread(file.path(p,'01_metadata/pretreatment_clinical_translation/ERP117672_PATIENT_RESPONSE_MAP.tsv'))[primary_included==TRUE,.(patient_id,group=ifelse(primary_group=='Responder','R','NR'))]
gsem <- fread(file.path(p,'01_metadata/phase15b_public_replication/GSE302495_PATIENT_RESPONSE_MAP.tsv'))[,.(patient_id,group=response_group)]
datasets <- c('ERP117672','GSE302495'); data <- list(); recovery<-list()
for(cohort in datasets){
  d<-fread(file.path(out,paste0('inputs/',cohort,'_MOLECULAR_INPUT_RESPONSE_BLIND.tsv')))
  mp<-if(cohort=='ERP117672') erpm else gsem
  d<-merge(d,mp,by='patient_id',all.x=TRUE,sort=FALSE)
  stopifnot(!anyNA(d),nrow(d)==if(cohort=='ERP117672')35 else 38,sum(d$group=='R')==if(cohort=='ERP117672')6 else 12)
  old<-fread(file.path(p,paste0('05_results/tables/phase15c_public_benchmark/',cohort,'_UNIVARIATE_BENCHMARK_EFFECTS.tsv')))
  for(v in c('PTPN22','AF35')){
    row<-old[variable==v]; g<-hedges_g(d[[v]],d$group); f<-fit_univariate_firth(d[[v]],d$group)
    loo<-vapply(seq_len(nrow(d)),function(i)hedges_g(d[[v]][-i],d$group[-i]),numeric(1))
    recovered<-c(g,f$or,f$low,f$high,min(loo),max(loo)); expected<-c(row$hedges_g,row$firth_or_per_sd,row$firth_ci_low,row$firth_ci_high,row$loo_hedges_g_min,row$loo_hedges_g_max)
    recovery[[paste(cohort,v)]]<-data.table(dataset=cohort,variable=v,metric=c('hedges_g','firth_or','firth_low','firth_high','loo_min','loo_max'),expected=expected,recovered=recovered,absolute_error=abs(recovered-expected),status=ifelse(abs(recovered-expected)<1e-8,'PASS','FAIL'))
  }
  data[[cohort]]<-d
  put(d,paste0('inputs/',cohort,'_PATIENT_ANALYSIS_INPUT.tsv'))
}
rec<-rbindlist(recovery);put(rec,'qa/FROZEN_IMPLEMENTATION_RECOVERY.tsv');stopifnot(all(rec$status=='PASS'))
writeLines(c('# 临床权威恢复审计','',
 'ANCHOR_READOUT_CLINICAL_AUTHORITY_GATE = PASS','',
 'Phase15C 完整清单及本轮几何锁哈希已核验。ERP 使用 40 人 log2(TPM+1) 空间逐基因 z 后的 AF35 等权均值，再正向 z；按冻结 eligibility 纳入 35 人（6 R /29 NR）。PTPN22 为相同 40 人空间表达 z。GSE 使用冻结 RDS molecular_scores 的 PTPN22_expression 和 AF35_score，纳入原始 38 人（12 R /26 NR），底层为 log2(CPM+1)。35 基因完整且不含 PTPN22；无重新选基因或校准。',
 '比较量直接来自 Phase15C RESPONSE_BLIND 表：broad_immune、CD274、IFNg_response、cytotoxicity、TCR_activation、checkpoint_exhaustion。程序为既有逐基因 z 的等权均值，CD274 为原始转换表达。相关系数不受正向仿射缩放影响。',
 '响应映射：ERP primary_included=TRUE 且 primary_group=Responder 编码 R，其余已冻结可评估者 NR；GSE response_group 原样 R/NR。以 patient_id 一对一连接，无患者池化或细胞级推断。',
 '直接提取并复用冻结 Phase15C hedges_g、safe_z、fit_univariate_firth 三个纯函数，不执行原工作流。4 个原始 g、12 个单变量 Firth OR/区间及 8 个 LOO 端点共 24 项恢复误差均 <1e-8。Bootstrap 为旧 20,000 次组内有放回方案的配对扩展，本轮两预测量严格共用索引；不要求与过去各自使用不同种子的边际 CI 数字一致。',
 'Firth 使用相同 logistf 默认优化控制、firth=TRUE 和 pl=TRUE。联合模型按本轮最终患者 universe 重新标准化两连续协变量。旧消融模块已有相同模型；后续核对属于重复呈现，不作独立证据。',
 '完整方法、固定种子和停止边界见 ANALYSIS_E_PROTOCOL_FREEZE.md；旧 checkpoint 的全局停止状态由用户本轮限定 Analysis E 授权扩展，不据此修改主稿或主图。'),file.path(out,'ANCHOR_READOUT_CLINICAL_AUTHORITY_AUDIT.md'),useBytes=TRUE)
summaries<-list(); allboot<-list(); allloo<-list(); joint<-list(); fits<-list()
for(k in seq_along(datasets)){
  cohort<-datasets[k]; d<-data[[cohort]]; n<-nrow(d); ir<-which(d$group=='R'); inn<-which(d$group=='NR')
  set.seed(2026090800L+k); draws<-20000L
  indices<-replicate(draws,c(sample(ir,length(ir),replace=TRUE),sample(inn,length(inn),replace=TRUE)))
  bg<-apply(indices,2,function(ii)c(PTPN22=hedges_g(d$PTPN22[ii],d$group[ii]),AF35=hedges_g(d$AF35[ii],d$group[ii])))
  delta<-bg[2,]-bg[1,]; stopifnot(all(is.finite(bg)))
  boot<-data.table(dataset=cohort,draw=seq_len(draws),g_PTPN22=bg[1,],g_AF35=bg[2,],delta_g=delta)
  put(data.table(draw=seq_len(draws),t(indices)),paste0('tables/',cohort,'_PAIRED_BOOTSTRAP_INDICES.tsv'))
  put(data.table(index=seq_len(n),patient_id=d$patient_id,group=d$group),paste0('tables/',cohort,'_INDEX_PATIENT_MAP.tsv'))
  loo<-rbindlist(lapply(seq_len(n),function(i)data.table(dataset=cohort,omitted_patient=d$patient_id[i],omitted_group=d$group[i],g_PTPN22=hedges_g(d$PTPN22[-i],d$group[-i]),g_AF35=hedges_g(d$AF35[-i],d$group[-i]))))
  loo[,delta_g:=g_AF35-g_PTPN22]
  gp<-hedges_g(d$PTPN22,d$group); ga<-hedges_g(d$AF35,d$group)
  summaries[[cohort]]<-data.table(dataset=cohort,n=n,n_R=length(ir),n_NR=length(inn),g_PTPN22=gp,g_AF35=ga,delta_g=ga-gp,delta_ci_low=quantile(delta,.025),delta_ci_high=quantile(delta,.975),bootstrap_positive_fraction=mean(delta>0),bootstrap_draws=draws,valid_draws=length(delta),seed=2026090800L+k,loo_delta_min=min(loo$delta_g),loo_delta_max=max(loo$delta_g),loo_positive_fraction=mean(loo$delta_g>0),g_PTPN22_ci_low=quantile(bg[1,],.025),g_PTPN22_ci_high=quantile(bg[1,],.975),g_AF35_ci_low=quantile(bg[2,],.025),g_AF35_ci_high=quantile(bg[2,],.975))
  allboot[[cohort]]<-boot;allloo[[cohort]]<-loo
  dd<-data.frame(y=as.integer(d$group=='R'),PTPN22_z=safe_z(d$PTPN22),AF35_z=safe_z(d$AF35))
  warn<-character(); fit<-withCallingHandlers(logistf(y~PTPN22_z+AF35_z,data=dd,firth=TRUE,pl=TRUE),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')})
  fits[[cohort]]<-fit
  convmax<-max(abs(fit$conv)); profmax<-max(fit$pl.iter); rankx<-qr(model.matrix(~PTPN22_z+AF35_z,dd))$rank
  bad<-!all(is.finite(c(coef(fit),fit$ci.lower,fit$ci.upper)))||convmax>=1e-3||profmax>=logistpl.control()$maxit||rankx<3
  for(v in c('PTPN22','AF35')){
    key<-paste0(v,'_z')
    joint[[paste(cohort,v)]]<-data.table(dataset=cohort,variable=v,n=n,n_R=sum(dd$y),n_NR=n-sum(dd$y),events_per_parameter=sum(dd$y)/2,beta=coef(fit)[key],OR=exp(coef(fit)[key]),beta_ci_low=fit$ci.lower[key],beta_ci_high=fit$ci.upper[key],OR_ci_low=exp(fit$ci.lower[key]),OR_ci_high=exp(fit$ci.upper[key]),profile_p=fit$prob[key],convergence_max=convmax,profile_iterations_max=profmax,profile_iterations_limit=logistpl.control()$maxit,design_rank=rankx,predictor_pearson=cor(dd$PTPN22_z,dd$AF35_z),predictor_VIF=1/(1-cor(dd$PTPN22_z,dd$AF35_z)^2),design_condition_number=kappa(model.matrix(~PTPN22_z+AF35_z,dd),exact=TRUE),warnings=if(length(warn))paste(unique(warn),collapse=';')else'NONE',gate=if(bad)'NOT_INTERPRETABLE'else if(sum(dd$y)/2<10)'PASS_WITH_LIMITS'else'PASS',interpretation='anchor/readout conditional-association sensitivity')
  }
  cat(cohort,'done\n')
}
put(rbindlist(summaries),'tables/AF35_VS_PTPN22_PAIRED_DELTA_G.tsv')
put(rbindlist(allboot),'tables/PAIRED_DELTA_G_BOOTSTRAP_DRAWS.tsv')
put(rbindlist(allloo),'tables/PAIRED_DELTA_G_PATIENT_LOO.tsv')
put(rbindlist(joint),'tables/PTPN22_AF35_JOINT_FIRTH_SENSITIVITY.tsv')
saveRDS(fits,file.path(out,'tables/JOINT_FIRTH_FITS.rds'))
capture.output(sessionInfo(),file=file.path(out,'logs/R_SESSION_INFO.txt'))
