options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table);library(digest)})
out <- 'AF35_MATCHED_RANDOM_CLINICAL_NULL_v1'; p <- 'HCC_ICI_project'
stopifnot(!file.exists(file.path(out,'AF35_MATCHED_NULL_PROGRAM_EFFECTS.tsv')))
locks <- rbind(fread(file.path(out,'MATCHED_PROGRAM_MANIFEST_SHA256.tsv')),fread(file.path(out,'SCORING_AUTHORITY_LOCK.tsv'))[,.(path,sha256)])
for(i in seq_len(nrow(locks))) stopifnot(digest(file.path(out,locks$path[i]),algo='sha256',file=TRUE)==locks$sha256[i])
manifest <- fread(file.path(out,'AF35_MATCHED_NULL_PROGRAM_MANIFEST.tsv'))
af <- fread(file.path(out,'AF35_ORDERED_GENES.tsv'))$gene
sets <- split(manifest$matched_gene,manifest$program_id)
stopifnot(length(sets)==2000L,all(lengths(sets)==35L),all(vapply(sets,uniqueN,integer(1))==35L))
if(!file.exists(file.path(out,'OUTCOME_JOIN_TIMESTAMP.txt'))) {
mats <- readRDS(file.path(out,'EXPRESSION_MATRICES_RESPONSE_BLIND.rds'))
z <- function(x) as.numeric((x-mean(x))/sd(x))
sc <- list()
for(cohort in names(mats)) {
  mat <- mats[[cohort]]; selected <- sort(unique(c(af,unlist(sets))))
  x <- mat[selected,,drop=FALSE]; stopifnot(all(is.finite(x)),all(apply(x,1,sd)>0))
  gz <- t(scale(t(x)))
  smat <- sapply(c(list(AF35=af),sets),function(g) z(colMeans(gz[g,,drop=FALSE])))
  stopifnot(all(is.finite(smat)),max(abs(colMeans(smat)))<1e-12,max(abs(apply(smat,2,sd)-1))<1e-12)
  raw_export <- data.table(gene=rownames(x),x)
  fwrite(raw_export,file.path(out,paste0(cohort,'_SELECTED_TRANSFORMED_EXPRESSION.tsv')),sep='\t')
  score_wide <- data.table(patient_id=colnames(mat),smat)
  sc[[cohort]] <- melt(score_wide,id.vars='patient_id',variable.name='program_id',value.name='score')[,cohort:=cohort]
}
scores <- rbindlist(sc)
fwrite(scores,file.path(out,'AF35_MATCHED_NULL_PATIENT_PROGRAM_SCORES_RESPONSE_BLIND.tsv'),sep='\t')
scorefiles <- c('AF35_MATCHED_NULL_PATIENT_PROGRAM_SCORES_RESPONSE_BLIND.tsv',paste0(names(mats),'_SELECTED_TRANSFORMED_EXPRESSION.tsv'))
fwrite(data.table(path=scorefiles,sha256=vapply(file.path(out,scorefiles),digest,character(1),algo='sha256',file=TRUE),time_utc=format(Sys.time(),tz='UTC',usetz=TRUE)),file.path(out,'PATIENT_SCORE_FREEZE_SHA256.tsv'),sep='\t')
# First outcome parsing occurs only below this persisted timestamp and hash checkpoint.
writeLines(format(Sys.time(),tz='UTC',format='%Y-%m-%dT%H:%M:%OS6Z'),file.path(out,'OUTCOME_JOIN_TIMESTAMP.txt'))
} else {
  score_lock <- fread(file.path(out,'PATIENT_SCORE_FREEZE_SHA256.tsv'))
  for(i in seq_len(nrow(score_lock))) stopifnot(digest(file.path(out,score_lock$path[i]),algo='sha256',file=TRUE)==score_lock$sha256[i])
  scores <- fread(file.path(out,'AF35_MATCHED_NULL_PATIENT_PROGRAM_SCORES_RESPONSE_BLIND.tsv'))
}
erp_path <- file.path(p,'01_metadata/pretreatment_clinical_translation/ERP117672_PATIENT_RESPONSE_MAP.tsv')
gse_path <- file.path(p,'01_metadata/phase15b_public_replication/GSE302495_PATIENT_RESPONSE_MAP.tsv')
stopifnot(digest(erp_path,algo='sha256',file=TRUE)=='4f28ae48fabf509f1248ac808d2c07b1a6beaf37b4f8f14fbeb999431f7082b2',digest(gse_path,algo='sha256',file=TRUE)=='fd2d8101719bd178e0b28680645e86319a30e485f210dd8f27fae04dc1c83b8d')
erp <- fread(erp_path); gse <- fread(gse_path)
stopifnot(!anyDuplicated(erp$patient_id),!anyDuplicated(gse$patient_id))
stopifnot(all(gse$response_group %in% c('R','NR')))
gse[,recist_group:=fifelse(response_group=='R','Responder','Nonresponder')]
maps <- rbind(erp[,.(cohort='ERP117672',patient_id,analysis_included=primary_included,response_group=primary_group)],gse[,.(cohort='GSE302495',patient_id,analysis_included=TRUE,response_group=recist_group)])
stopifnot(maps[cohort=='ERP117672' & analysis_included==TRUE,.N]==35L,maps[cohort=='ERP117672' & analysis_included==TRUE & response_group=='Responder',.N]==6L,maps[cohort=='GSE302495' & response_group=='Responder',.N]==12L,nrow(maps)==78L)
joined <- merge(scores,maps,by=c('cohort','patient_id'),all.x=TRUE,sort=FALSE)
stopifnot(nrow(joined)==2001L*78L,!anyNA(joined$analysis_included))
fwrite(joined,file.path(out,'AF35_MATCHED_NULL_PATIENT_PROGRAM_SCORES.tsv'),sep='\t')
fwrite(maps,file.path(out,'ANALYSIS_PATIENT_ELIGIBILITY.tsv'),sep='\t')
hg <- function(x,gr) {
  r <- x[gr=='Responder']; nr <- x[gr=='Nonresponder']; n <- length(r)+length(nr)
  s <- sqrt(((length(r)-1)*var(r)+(length(nr)-1)*var(nr))/(n-2))
  (1-3/(4*(n-2)-1))*(mean(r)-mean(nr))/s
}
effects <- joined[analysis_included==TRUE,.(n=.N,n_R=sum(response_group=='Responder'),n_NR=sum(response_group=='Nonresponder'),hedges_g=hg(score,response_group)),by=.(cohort,program_id)]
stopifnot(nrow(effects)==4002L,all(is.finite(effects$hedges_g)))
effects[,direction:=fifelse(hedges_g>0,'POSITIVE',fifelse(hedges_g<0,'NEGATIVE','ZERO'))]
frozen_erp <- fread(file.path(p,'05_results/tables/pretreatment_clinical_translation/ERP117672_AF_PROGRAM_RESPONSE_ASSOCIATION.tsv'))[population=='PRIMARY_RESPONSE_EVALUABLE',hedges_g]
frozen_gse <- fread(file.path(p,'05_results/tables/phase15b_public_replication/GSE302495_AF35_RECIST_ASSOCIATION.tsv'))$hedges_g
reproduction <- effects[program_id=='AF35']
reproduction[,frozen_hedges_g:=ifelse(cohort=='ERP117672',frozen_erp,frozen_gse)]
reproduction[,absolute_error:=abs(hedges_g-frozen_hedges_g)]
reproduction[,pass:=absolute_error<1e-10]
fwrite(reproduction,file.path(out,'AF35_FULL_SCORE_EFFECT_REPRODUCTION.tsv'),sep='\t')
stopifnot(all(reproduction$pass))
fwrite(effects,file.path(out,'AF35_MATCHED_NULL_PROGRAM_EFFECTS.tsv'),sep='\t')
cat('AF35_MATCHED_NULL_EFFECT_GATE=PASS\n');print(reproduction)
