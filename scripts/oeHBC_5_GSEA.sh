#!/bin/bash 


ncores=23 
nrot=10000000
NOW=$(date +"%m%d%Y_%H%M%S")
module load R-3.3.1

R --vanilla <<< "rmarkdown::render('oeHBCregenWT_GSEAprep.Rmd')"

R --vanilla < oeHBCregenWT_romerGSEA.R --args --expt oeHBCregenWT --ncores $ncores --nrot $nrot > "oeHBCregenWT_GSEA_"$NOW.Rout
