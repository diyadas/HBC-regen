#title: Make SummarizedExperiment objects from SCONE output
#author: Diya Das and Russell Fletcher
#date: November 11, 2017

rm(list=ls()); options(getClass.msg=FALSE)
library(scone)
library(BiocParallel)
library(optparse)
library(SummarizedExperiment)


option_list <- list(
  make_option("--expt", default="", type="character", help="full form, e.g. Expt1"),
  make_option("--ncores", default="1", type="double"),
  make_option("--norm", type="character")
)

opt <- parse_args(OptionParser(option_list=option_list))
expt_str <- opt$expt

register(MulticoreParam(workers = opt$ncores))

out_dir <- paste0("../output/clust/",expt_str)
scone_out <- get(load(file.path(out_dir, paste0(expt_str,"_", opt$norm, "_scone.Rda")))[1])

load(file.path(out_dir,paste0(expt_str,"_filtdata.Rda")))

norm <- opt$norm
normstr <- gsub(",","_",gsub("_k=|_","",norm))
scone_out <- scone_out$normalized_data[[norm]]
se <- SummarizedExperiment(list(counts = expm1(scone_out)),
                           colData=data.frame(batch=batch, expt=expt,  nreads=qc$NREADS, ralign=qc$RALIGN, pctmrna=qc$PCT_MRNA_BASES))
save(se, file = file.path(out_dir,paste0(expt_str,"_",normstr,"_se.Rda")))
