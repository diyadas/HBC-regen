#title: Excluding Contaminant Cells
#author: Diya Das and Russell Fletcher
#date: November 11, 2017

library(SummarizedExperiment)
out_dir <- "../output/clust/oeHBC"

normfiles=list.files(path = out_dir, pattern=glob2rx("oeHBC_*fq*_se.Rda"))

reg3g_list_all = list()
omp_cyp_list_all = list()
csf1r_eGFP_list_all = list()

lapply(normfiles, function(sefile){
norm <- limma::removeExt(sefile) 

load(file.path(out_dir, sefile))
counts <- assay(se)
expt <- colData(se)$expt

omp_cyp <- counts["Omp",] > 25 & (counts["Cyp2g1",] > 25 | counts["Cyp1a2",] > 25)
omp_cyp_list <- colnames(counts)[omp_cyp]

csf1r_eGFP <- counts["Csf1r",] > 25 & counts["eGFP",] < 25
csf1r_eGFP_list <- colnames(counts)[csf1r_eGFP]

reg3g_list <- colnames(counts[,grep("UI|SOX2EGFP",expt)])[counts["Reg3g",grep("UI|SOX2EGFP",expt)]>=25]

save(reg3g_list,omp_cyp_list, csf1r_eGFP_list,file=file.path("../ref",paste0(norm,"_exclude.Rda")))

#rm(list=setdiff(ls(), c('out_dir','normfiles')))
reg3g_list_all <<- union(reg3g_list_all,reg3g_list)
omp_cyp_list_all <<- union(omp_cyp_list_all, omp_cyp_list)
csf1r_eGFP_list_all <<- union(csf1r_eGFP_list_all, csf1r_eGFP_list)
})
save(reg3g_list_all, omp_cyp_list_all, csf1r_eGFP_list_all, file=file.path("../ref","oeHBC_combined_exclude.Rda"))