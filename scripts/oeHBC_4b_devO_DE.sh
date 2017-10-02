#!/bin/bash  

NOW=$(date +"%m%d%Y_%H%M%S")

module load R-3.3.1

R --vanilla<<<"rmarkdown::render('oeHBCregenWT_slingshot.Rmd')"
R --vanilla<<<"rmarkdown::render('oeHBCregenWTKO_slingshot.Rmd')"
R --vanilla<<<"rmarkdown::render('oeHBCregenWT_DE.Rmd')"
R --vanilla<<<"rmarkdown::render('oeHBCregenWTKO_DE.Rmd')"