#!/bin/bash 
#
###Download and install homer.
###Place list of the top 1000 genes enriched in the HBC*1 cluster relative to the resting HBCs (HBC) cluster
###in the homer directory ("HBC1top1000LFC.txt")

findMotifs.pl HBC1_top1000LFC.txt mouse /output/DE/oeHBCregenWT -start -1000 -end 100 -len 8,16 -p 4