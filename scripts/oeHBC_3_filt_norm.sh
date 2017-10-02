#!/bin/bash 


ncores=11
 
NOW=$(date +"%m%d%Y_%H%M%S")
norm="none,fq,ruv_k=1,no_bio,no_batch"
module load R-3.3.1

R --vanilla < oeHBC_filtering.R --args --expt oeHBC --exclude combined > "oeHBC_filtering_"$NOW.Rout
R --vanilla < oeHBC_norm.R --args --expt oeHBC --ncores $ncores --norm $norm  > "oeHBC_norm_"$NOW.Rout
R --vanilla < oeHBC_makeSE_expt.R --args --expt oeHBC --ncores $ncores --norm $norm  > "oeHBC_makeSE_expt_"$NOW.Rout
