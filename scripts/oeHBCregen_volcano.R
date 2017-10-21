# title: Volcano Plots for Lineage Transitions in Regeneration
# author: Diya Das
# date: April 6, 2017


library(clusterExperiment); library(calibrate); library(plyr)
expt_str <- "oeHBCregenWT"
DE_dir <- file.path("../output/DE", expt_str)
clust_dir <- file.path("../output/clust", expt_str)

tfs <- unlist(read.table("../ref/ATFDB_mm_TF.txt"))

load(file.path(clust_dir, paste0(expt_str,"_cmmerged.Rda")))
load(file.path(clust_dir, paste0(expt_str, "_slingshot_out.Rda")))

pairsDE <- getBestFeatures(transform(cmobj), primaryCluster(cmobj), contrastType="Pairs",number=Inf,p.value=0.05)
write.table(pairsDE, file.path(DE_dir,"oeHBCregenWT_pairwiseDE.txt"),quote=FALSE,sep="\t")

get_tuple <- function(lineages, n) rbind(do.call(rbind,lapply(lineages[n], function(x) data.frame(cbind(x[-1],x[-length(x)]), stringsAsFactors = F))))

contrast_tuple <- unique(get_tuple(lineages, 1:4))
contrast_tuple <- rbind(resting1=c(9,1), contrast_tuple)
contrast_tuple[contrast_tuple==11] <- 10
contrast_tuple[contrast_tuple==12] <- 11

contrast_name <- data.frame(contrast = paste0("X",contrast_tuple$X1,"-X",contrast_tuple$X2), names= c("H1-H", "H2-H1", "G-H2", "I12-G", "I3-I12", "iO-I3", "mO-iO", "R1-H2","H-R1", "R2-R1", "S-H2"))

genes_to_label = c("Trp63","Krt5","Krt14","Cyp2g1","Notch2","Cyp1a2","Trpm5","Ascl3","Cftr","Sox9", "Kit","Ascl1","Ccnd1","Top2a","Neurod1","Lhx2","Gap43","Gng8","Omp","Cnga2","Gng13", "Reg3g", "Krtdap","Sprr1a", "Krt16")

label_markers <- function(cont_sub, pval_thresh, genes_to_label){
  # function to label top 6 genes + selected markers
  with(subset(cont_sub, adj.P.Val<pval_thresh & abs(logFC)>1), {
    n = unique(c(na.omit(order(logFC)[1:6]),
                 na.omit(order(logFC, decreasing = T)[1:6]),
                 na.omit(match(genes_to_label,Feature))))
    n1 <- vector(); n2 <- vector()
    invisible(lapply(n, function(z) {
      if (logFC[z] > 0) {
        n1 <<- c(n1,z)
      }else{
        n2 <<- c(n2,z)
      }
    }))
    points(logFC[n1], -log10(adj.P.Val[n1]), pch=1, col=c("firebrick4"), cex = 0.7)
    textxy(logFC[n1], -log10(adj.P.Val[n1]), labs=Feature[n1], cex=.5, offset=0.7, col=c("firebrick4"))
    points(logFC[n2], -log10(adj.P.Val[n2]), pch=1, col=c("dodgerblue4"), cex = 0.7)
    textxy(logFC[n2], -log10(adj.P.Val[n2]), labs=Feature[n2], cex=.5, offset=0.7, col=c("dodgerblue4"))
  })
}

tmplist <- apply(contrast_tuple, 1, function(x){
  desired_contrast <- paste0("X",x[1],"-X",x[2])
  opposite_contrast <- paste0("X",x[2],"-X",x[1])
  return(ifelse (desired_contrast %in% levels(pairsDE$Contrast), desired_contrast, opposite_contrast))})
pairsDE_sub <- pairsDE[pairsDE$Contrast %in% tmplist,]

xlim_1 <- round_any(max(abs(pairsDE$logFC)), 1, f=ceiling)
ylim_2 <- round_any(max(-log10(pairsDE_sub$adj.P.Val)), 10, f = ceiling)

plotContrasts <- function(DE, contrast_tuple, contrast_name, fname, pval_thresh){
  write(c("Contrast name", "Positive Change", "Negative change", paste("max_pval <", pval_thresh)), file=fname, sep = "\t", append=FALSE, ncolumns = 4)
  tmplist <- apply(contrast_tuple, 1, function(x){
    desired_contrast <- paste0("X",x[1],"-X",x[2])
    opposite_contrast <- paste0("X",x[2],"-X",x[1])
    if (desired_contrast %in% levels(DE$Contrast)){
      cont_sub <<- subset(DE, abs(logFC) > 1 & Contrast==desired_contrast)
    } else if (opposite_contrast %in% levels(DE$Contrast)) {
      cont_sub <<- subset(DE, abs(logFC) > 1 & Contrast==opposite_contrast)
      cont_sub$logFC <<- -cont_sub$logFC
    } else {
      stop('stop something wrong')
    }
    
    # plotting all genes
    png(file.path(DE_dir,paste0(expt_str, "_", contrast_name$names[contrast_name$contrast==desired_contrast],"_volcano_",Sys.Date(),".png")), width = 2.25, height = 2.25, units="in", res=600)
    par(family='Helvetica', cex.lab=0.5, cex.axis=0.5, oma=c(0,0,0,0), mgp=c(0.8,0.15,0), xpd=F, lwd=0.5, pty="s", mar=c(1.7,1.7,0.3,0.3))
    with(cont_sub, plot(logFC, -log10(adj.P.Val),xlim=c(-xlim_1,xlim_1), pch=20, col="grey", ylim=c(0,90),cex=0.5, xaxt='n',yaxt='n',xlab='logFC',ylab='-log10(adj. p-value)'))
    axis(1, lwd = 0.5, tck=-0.01*2, labels=T);axis(2, lwd = 0.5, tck=-0.01*2, labels=T)
    with(subset(cont_sub, adj.P.Val<pval_thresh & abs(logFC)>1), points(logFC, -log10(adj.P.Val), pch=20, col="black", cex=0.5))
    
    # label number of DE genes
    neg_num <- with(subset(cont_sub, adj.P.Val<pval_thresh & logFC < -1), length(Feature))
    pos_num <- with(subset(cont_sub, adj.P.Val<pval_thresh & logFC > 1), length(Feature))
    usr <- par( "usr" )
    text( usr[ 1 ], usr[ 4 ], neg_num, adj = c( -0.5, 2 ), col = "dodgerblue4", font=4)
    text( usr[ 2 ], usr[ 4 ], pos_num, adj = c( 1.5, 2 ), col = "firebrick4", font=4)
    
    # label_markers(cont_sub, pval_thresh, genes_to_label) # sanity check
    dev.off()
    
    # plotting only transcription factors
    png(file.path(DE_dir,paste0(expt_str, "_", contrast_name$names[contrast_name$contrast==desired_contrast],"_tf_volcano_",Sys.Date(),".png")), width = 2.25, height = 2.25, units="in", res=600)
    par(family='Helvetica', cex.lab=0.5, cex.axis=0.5, oma=c(0,0,0,0), mgp=c(0.8,0.15,0),xpd=F, lwd=0.5, pty="s", mar=c(1.7,1.7,0.3,0.3))
    
    with(subset(cont_sub, Feature %in% tfs), plot(logFC, -log10(adj.P.Val),xlim=c(-xlim_1,xlim_1), pch=20, col="grey", ylim=c(0,90),cex=0.5, xaxt='n',yaxt='n',xlab='logFC',ylab='-log10(adj. p-value)'))
    axis(1, lwd = 0.5, tck=-0.01*2, labels=T);axis(2, lwd = 0.5, tck=-0.01*2, labels=T)
    with(subset(cont_sub, adj.P.Val<pval_thresh & abs(logFC)>1 & Feature %in% tfs), points(logFC, -log10(adj.P.Val), pch=20, col="black", cex=0.5))
    
    #label number of DE genes
    neg_num <- with(subset(cont_sub, adj.P.Val<pval_thresh & logFC < -1 & Feature %in% tfs), length(Feature))
    pos_num <- with(subset(cont_sub, adj.P.Val<pval_thresh & logFC > 1 & Feature %in% tfs), length(Feature))
    usr <- par( "usr" )
    text( usr[ 1 ], usr[ 4 ], neg_num, adj = c( -0.5, 2 ), col = "dodgerblue4", font=4)
    text( usr[ 2 ], usr[ 4 ], pos_num, adj = c( 1.5, 2 ), col = "firebrick4", font=4)
    
    dev.off()
    
    with(subset(cont_sub, adj.P.Val<pval_thresh & abs(logFC)>1), write(c(contrast_name$names[contrast_name$contrast==desired_contrast], sum(logFC>0), sum(logFC<0), max(adj.P.Val)), ncolumns=4, file=fname, sep = "\t", append=TRUE))
  })
}  

pval_thresh <- 1e-2

plotContrasts(pairsDE, contrast_tuple, contrast_name, fname = file.path(DE_dir,paste0(expt_str,"_pairsDE_tuple_p",pval_thresh,"_FC1_2_",Sys.Date(),".txt")), pval_thresh)
