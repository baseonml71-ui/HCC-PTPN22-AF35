suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
root <- 'AF35_LONGITUDINAL_CLONAL_STATE_DYNAMICS_M2_CLEANROOM_v1'
f <- file.path(root,'02_phaseA_frozen_outputs'); b <- file.path(root,'03_phaseB_outcome_inputs'); cdir <- file.path(root,'04_phaseB_results')
rd <- function(dir,name) fread(file.path(dir,name),header=TRUE)
checks <- list()
check <- function(name,a,b,tol=1e-10) {
  ok <- length(a)==length(b) && identical(as.vector(is.na(a)),as.vector(is.na(b)))
  err <- if(ok && is.numeric(a) && is.numeric(b) && any(!is.na(a))) max(abs(a-b),na.rm=TRUE) else if(ok && identical(as.character(a),as.character(b))) 0 else Inf
  checks[[length(checks)+1L]] <<- data.table(check=name,max_abs_error=err,pass=is.finite(err)&&err<=tol)
}
lock <- fromJSON(file.path(root,'AF35_M2_PHASEA_FREEZE_LOCK.json'))
for(n in names(lock$hashes)) check(paste('frozen hash',n),digest(file=file.path(root,n),algo='sha256'),lock$hashes[[n]])
check('phaseA manifest hash',digest(file=file.path(root,'AF35_M2_PHASEA_COMPLETE_SHA256_MANIFEST.tsv'),algo='sha256'),lock$manifest_sha256)
start <- fromJSON(file.path(root,'QA/PHASEB_START_EVENT.json'))
check('join starts after freeze',start$utc>lock$utc,TRUE)
lab <- rd(b,'RESPONSE.tsv');check('clinical input fields',names(lab),c('patient','response'));check('no duplicate labels',uniqueN(lab$patient),nrow(lab))
delta <- rd(f,'AF35_M2_PATIENT_LONGITUDINAL_CHANGES_RESPONSE_BLIND.tsv'); j <- rd(cdir,'PATIENT_CHANGES_AFTER_JOIN.tsv')
check('joined patient set',sort(j$patient),sort(delta$patient));check('joined labels',j$response,lab$response[match(j$patient,lab$patient)])
for(col in names(delta)) check(paste('joined unchanged',col),j[[col]],delta[[col]][match(j$patient,delta$patient)])
verify <- function(prefix,values,cols,names0) {
  ord <- rd(cdir,paste0(prefix,'_PATIENT_ORDER.tsv')); ii <- match(ord$patient,values$patient); values <- values[ii]
  v <- as.matrix(values[,cols,with=FALSE]);ir <- which(ord$response=='R');inn <- which(ord$response=='NR')
  obs <- colMeans(v[ir,,drop=FALSE])-colMeans(v[inn,,drop=FALSE])
  idx <- as.matrix(rd(cdir,paste0(prefix,'_BOOTSTRAP_INDICES.tsv')))+1L
  check(paste(prefix,'10000 draws'),nrow(idx),10000)
  check(paste(prefix,'R stratification'),all(idx[,seq_along(ir),drop=FALSE] %in% ir),TRUE)
  check(paste(prefix,'NR stratification'),all(idx[,length(ir)+seq_along(inn),drop=FALSE] %in% inn),TRUE)
  sims <- t(apply(idx,1,function(i) colMeans(v[i[seq_along(ir)],,drop=FALSE])-colMeans(v[i[length(ir)+seq_along(inn)],,drop=FALSE])))
  if(length(cols)==1) sims <- matrix(sims,ncol=1)
  use <- names0
  if(length(cols)==3) {sims <- cbind(sims,sims[,2]-sims[,3]);obs <- c(obs,obs[2]-obs[3]);use <- c(use,'THETA_BALANCE')}
  boot <- rd(cdir,paste0(prefix,'_BOOTSTRAP.tsv'));summary <- rd(cdir,paste0(prefix,'_BOOTSTRAP_SUMMARY.tsv'))
  for(k in seq_along(use)) {
    row <- summary[component==use[k]]
    check(paste(prefix,'estimate',use[k]),obs[k],row$estimate)
    check(paste(prefix,'all bootstrap',use[k]),sims[,k],boot[[use[k]]])
    check(paste(prefix,'CI',use[k]),as.vector(quantile(sims[,k],c(.025,.975),type=7)),c(row$bootstrap_ci_low,row$bootstrap_ci_high))
  }
  assignments <- combn(nrow(v),length(ir));dist <- t(apply(assignments,2,function(i) colMeans(v[i,,drop=FALSE])-colMeans(v[-i,,drop=FALSE])))
  if(length(cols)==1) dist <- matrix(dist,ncol=1)
  py <- rd(cdir,paste0(prefix,'_EXACT_PERMUTATION_DISTRIBUTION.tsv'))
  rkeys <- apply(assignments-1L,2,paste,collapse=',');matchi <- match(py$R_indices,rkeys)
  check(paste(prefix,'complete exact assignments'),sort(py$R_indices),sort(rkeys))
  for(k in seq_along(names0)) {
    check(paste(prefix,'exact distribution',names0[k]),dist[matchi,k],py[[names0[k]]])
    check(paste(prefix,'exact p',names0[k]),mean(abs(dist[,k])>=abs(obs[k])-1e-12),summary[component==names0[k],exact_two_sided_p])
  }
  if(length(cols)==3) {check('theta additive identity',obs[1],obs[2]+obs[3]);check('bootstrap additive identity',sims[,1],sims[,2]+sims[,3])}
}
verify('AF35_M2',j,c('CHANGE_O_CELL','CHANGE_O_CLONE','CHANGE_D_SIZE'),c('THETA_CELL','THETA_CLONE','THETA_SIZE'))
shared <- rd(f,'AF35_M2_SHARED_CLONE_METRICS_RESPONSE_BLIND.tsv')[shift_eligible==TRUE]
if(file.exists(file.path(cdir,'AF35_M2_SHARED_CLONE_CLINICAL.tsv'))) verify('AF35_M2_SHARED',shared,'SHARED_CLONE_MEAN_SHIFT','THETA_SHARED')
sc <- rd(f,'AF35_M2_STATE_CONTRIBUTIONS_RESPONSE_BLIND.tsv');sc <- merge(sc,lab,by='patient',all.x=TRUE)
out <- rd(cdir,'STATE_CONTRIBUTIONS.tsv')
for(gr in c('ALL','R','NR','R_MINUS_NR')) for(st in unique(out$released_state)) {
  z <- sc[released_state==st];py <- out[group==gr & released_state==st]
  for(col in c('CELL_CONTRIBUTION','CLONE_CONTRIBUTION')) {
    value <- if(gr=='ALL') mean(z[[col]]) else if(gr=='R_MINUS_NR') mean(z[response=='R'][[col]])-mean(z[response=='NR'][[col]]) else mean(z[response==gr][[col]])
    check(paste('state',gr,st,col),value,py[[col]])
  }
}
comp <- rd(f,'AF35_M2_PATIENT_TIMEPOINT_COMPONENTS_RESPONSE_BLIND.tsv');comp <- merge(comp,lab,by='patient',all.x=TRUE)
desc <- rd(cdir,'WITHIN_GROUP_LONGITUDINAL_DESCRIPTION.tsv')
for(i in seq_len(nrow(desc))) {
  row <- desc[i];p <- row$response;col <- row$component;v <- j[response==p][[paste0('CHANGE_',col)]]
  check(paste('description',p,col),c(mean(comp[response==p & timepoint=='pre'][[col]]),mean(comp[response==p & timepoint=='post'][[col]]),mean(v),median(v),sum(v>0)),as.numeric(row[,.(pre_mean,later_mean,change_mean,change_median,positive_change_patients)]))
}
# Baseline is compared with its frozen numeric source, not re-estimated.
source <- rd('AF35_CLONOTYPE_NEUTRALIZED_STATE_OCCUPANCY_v1','AF35_CLONE_NEUTRALIZATION_CLINICAL_DECOMPOSITION.tsv')
base <- rd(b,'L_FROZEN_GROUP_BASELINE.tsv')
for(col in names(base)) check(paste('frozen baseline copy',col),base[[col]],source[[col]][match(base$component,source$component)])
res <- rbindlist(checks);fwrite(res,file.path(root,'QA/INDEPENDENT_R_PHASEB_CHECKS.tsv'),sep='\t');writeLines(capture.output(sessionInfo()),file.path(root,'QA/R_PHASEB_SESSION_INFO.txt'))
print(res[pass==FALSE]);stopifnot(all(res$pass));cat(nrow(res),'Phase B independent checks passed; max error:',max(res$max_abs_error),'\n')
