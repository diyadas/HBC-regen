#!/bin/bash  

NOW=$(date +"%m%d%Y_%H%M%S")

module load R-3.3.1

R --vanilla <<< "rmarkdown::render('oeHBCregen_clusterPlots.Rmd')"
R --vanilla <<< "rmarkdown::render('oeHBCdiffregen_clusterPlots.Rmd')"
R --vanilla <<< "rmarkdown::render('oeHBCregen_devorderplots.Rmd')"
R --vanilla <<< "rmarkdown::render('oeHBCregen_tf.Rmd')"
R --vanilla <<< "rmarkdown::render('oeHBCregenWT_genePairsPlots.Rmd')"
R --vanilla <<< "rmarkdown::render('oeHBCregenWT_genePlots.Rmd')"
R --vanilla <<< "rmarkdown::render('oeHBCregenWT_AP1_TF_Heatmaps.Rmd')"
R --vanilla < oeHBCregen_volcano.R
R --vanilla < oeHBCregen_OR.R
