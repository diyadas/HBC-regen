#title: Make SummarizedExperiment objects from SCONE output
#author: Diya Das and Russell Fletcher
#date: November 11, 2017

rm(list=ls()); options(getClass.msg=FALSE)
library(scone)
library(BiocParallel)
library(optparse)
library(SummarizedExperiment)

####----load post filtering Rda to get expt, batch, QC
load(paste0("../output/clust/oeHBC/oeHBC_filtdata.Rda"))
names(batch) <- names(expt) <-rownames(qc)

####----load scone output (normalized data matrix) Rda
load(paste0("../output/clust/oeHBC/oeHBC_scone.Rda"))

expt_list <- c("oeHBCregenWT","oeHBCregenWTKO","oeHBCdiffregen")

lapply(names(scone_out$normalized_data), function(nrm){
nrmstr = gsub(",","_",gsub("_k=|_","",nrm))

counts_all <- scone_out$normalized_data[[nrm]]

lapply(expt_list, function (expt_str) {
if (expt_str == "oeHBCregenWTKO") {
  desiredSamples = gdata::startsWith(expt, "K5ERRY_", ignore.case=T) | gdata::startsWith(expt, "K5ERSox2cKO_", ignore.case=T)
}else if(expt_str == "oeHBCregenWT") {
  desiredSamples = gdata::startsWith(expt, "K5ERRY_", ignore.case=T)
}else if(expt_str == "oeHBCdiffregen") {
  desiredSamples = gdata::startsWith(expt, "K5ERRY_", ignore.case=T) | gdata::startsWith(expt, "K5ERP63CKO_", ignore.case=T)  | gdata::startsWith(expt, "SOX2EGFP+ICAM-F3-SRB1-", ignore.case=T)
}else {
  stop("code should never be reached")
}

counts <- counts_all[, desiredSamples]
out_dir <- paste0("../output/clust/",expt_str)
esh <- paste0("E",substr(expt_str,5,nchar(expt_str)))
se <- SummarizedExperiment(list(counts = expm1(counts)),
                           colData=data.frame(batch=batch[colnames(counts)], expt=expt[colnames(counts)],  nreads=qc[colnames(counts),]$NREADS, ralign=qc[colnames(counts),]$RALIGN, pctmrna=qc[colnames(counts),]$PCT_MRNA_BASES))
save(se, file = file.path(out_dir,paste0(esh,"_",nrmstr,"_se.Rda")))
})

})