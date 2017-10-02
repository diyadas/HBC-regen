#!/bin/bash 


ncores=23
 
NOW=$(date +"%m%d%Y_%H%M%S")
norm="none,fq,ruv_k=1,no_bio,no_batch"
norm2="none,fq,ruv_k=1,bio,no_batch"
norm3="none,fq,ruv_k=1,bio,batch"

module load R-3.3.1


R --vanilla < oeHBC_norm.R --args --expt oeHBC --ncores $ncores --norm $norm
R --vanilla < oeHBC_makeSE.R --args --expt oeHBC --ncores $ncores --norm $norm
R --vanilla < oeHBC_norm.R --args --expt oeHBC --ncores $ncores --norm $norm2
R --vanilla < oeHBC_makeSE.R --args --expt oeHBC --ncores $ncores --norm $norm2
R --vanilla < oeHBC_norm.R --args --expt oeHBC --ncores $ncores --norm $norm3
R --vanilla < oeHBC_makeSE.R --args --expt oeHBC --ncores $ncores --norm $norm3
R --vanilla < oeHBC_exclude.R
