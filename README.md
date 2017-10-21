# HBC-regen
Code and resources related to the analysis of regeneration in the olfactory epithelium

Below are the R scripts for analyzing the single-cell RNA-seq data from HBC stem cells of the olfactory epithelium, presented in the following manuscript:

Gadye L\*, Das D\*, Sanchez MA\*, Street KN, Baudhuin A, Wagner A, Cole MB, Choi YG, Yosef N, Purdom E, Dudoit S, Risso D, Ngai J, Fletcher RB. Identification of activated stem cell states unique to regeneration in the olfactory epithelium. (\* equal contribution)

The data are available on GEO in [GSE99251](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99251) and [GSE95601](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95601).

The repository currently has scripts that take as input Expression Set data and perform a series of computations, interspersed with visualizations. First, the data are filtered for poor quality cells and less informative genes. The data are normalized, and biological contaminants and known doublets (based on co-expression of differentiated cell markers) are removed. Then, the data are re-filtered and re-normalized.

After filtering and normalization, we clustered the data using clusterExperiment, performed developmental ordering and inferred lineage trajectories and branching with slingshot. For each lineage, differentially expressed genes were identified. We used Gene Set Enrichment Analysis to infer pathways regulating cell fates and transitions.

We created a number of visualizations based on clustering, experimental condition, and developmental order. We displayed coordinated and correlated differentially expressed genes including transcription factors, as well as a set of cell cycle genes and selected regulators of cell fate transitions along each lineage. The olfactory receptors and factors associated with OR regulation were plotted along the neuronal lineage. We also presented the top enriched gene sets for each cell cluster. 


### Create output directories and add to .gitignore
In project directory, run `mkdir -p output/{clust,data,romer,viz,DE,EDA}/oeHBCregen`, and add new directories to `.gitignore`. Place the scripts in the 'scripts' directory and the initial eSet 'data' in the data directory.

### Filtering and Normalization
`oeHBC_1_filt_norm.sh` performs the following analyses, by calling various R scripts (given in parentheses):

1. Filtering based on technical attributes (`oeHBC_filtering.R`)
2. Normalization using SCONE to get rankings (`oeHBC_norm.R`)

`oeHBC_2_norm.sh` performs the following analyses, by calling various R scripts (given in parentheses):

3. Get normalized matrices for several normalizations (`oeHBC_norm.R`)
3. Make SummarizedExperiment objects for each normalization (`oeHBC_makeSE.R`)
4. Create final list of samples to exclude as biological contaminants (`oeHBC_exclude.R`)

`oeHBC_3_filt_norm.sh` performs the following analyses, by calling various R scripts:

5. Filtering based on contaminants and technical attributes (`oeHBC_filtering.R`)
4. Re-normalization after removal of contaminants (`oeHBCregen_norm.R`)
5. Create SummarizedExperiment object for each experiment (`oeHBCregen_makeSE_expt.R`)



### Clustering, Developmental Ordering, Differential Expression
`oeHBC_4_clust.sh` performs the following:

6. Cluster samples in each experiment (`oeHBC_clust.R`)

`oeHBC_4b_devO_DE.sh` performs the following analyses, by calling various R scripts (given in parentheses):

8. Developmental ordering with slingshot (`oeHBCregenWT_slingshot.Rmd` & `oeHBCregenWTKO_slingshot.Rmd`)
9. Differential gene expression using limma, along each lineage (`oeHBCregenWT_DE.Rmd`)

`oeHBC_5_GSEA.sh`:

11. Preparation of gene sets for Gene Set Enrichment Analysis (GSEA; `oeHBCregenWT_GSEAprep.Rmd`)
11. GSEA based on cell clustering using limma romer (`oeHBCregenWT_romerGSEA.R`)

### Visualizations
`oeHBC_6_viz.sh` performs the following analyses, by calling various R scripts (given in parentheses):

13. Visualizations based on cell clustering (heatmap of marker genes, tSNE plots, PCA pairs plot, cluster & experimental condition bubble plots; `oeHBCregen_clusterPlots.Rmd` & `oeHBCdiffregen_clusterPlots.Rmd`)
12. Visualizations incorporating developmental ordering (3D-PCA plots, dot plots; `oeHBCregen_devorderplots.Rmd`)
12. Transcription factor co-expression, network analysis, and visualizations (`oeHBCregen_tf.Rmd`)
13. Plots of individual or pairs of genes in developmental order (`oeHBCregen_genePlots.Rmd` & `oeHBCregen_genePairsPlots.Rmd`)
15. Volcano plots of differentially expressed genes (`oeHBCregen_volcano.R`)
16. Olfactory Receptor (OR) gene and OR regulation associated gene expression plots (`oeHBCregen_OR.R`)

### Motif Analysis
17. Look for transcription factor motifs in the top 1000 most enriched genes in the activated HBC1 cluster relative to resting HBCs (`oeHBC1_findMotif.sh`)

### Dependencies/useful R packages:

- SCONE (normalization): https://github.com/YosefLab/scone
- clusterExperiment (clustering): https://github.com/epurdom/clusterExperiment
- slingshot (lineage trajectory algorithm): https://github.com/kstreet13/slingshot
