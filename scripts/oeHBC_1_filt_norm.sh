#!/bin/bash 


ncores=11
 
NOW=$(date +"%m%d%Y_%H%M%S")
module load R-3.3.1

R --vanilla < oeHBC_filtering.R --args --expt oeHBC > "oeHBC_filtering_"$NOW.Rout
R --vanilla < oeHBC_norm.R --args --expt oeHBC --ncores $ncores  > "oeHBC_norm_"$NOW.Rout
