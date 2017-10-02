#!/bin/bash 


ncores=23 
NOW=$(date +"%m%d%Y_%H%M%S")
norm="none,fq,ruv_k=1,no_bio,no_batch"
module load R-3.3.1

R --vanilla < oeHBC_clust.R --args --expt oeHBCregenWT --ncores $ncores --norm $norm > "oeHBC_clust_"$NOW.Rout
R --vanilla < oeHBC_clust.R --args --expt oeHBCregenWTKO --ncores $ncores --norm $norm > "oeHBC_clust_"$NOW.Rout
R --vanilla < oeHBC_clust.R --args --expt oeHBCdiffregen --ncores $ncores --norm $norm > "oeHBC_clust_"$NOW.Rout
