exfiles = list.files(path = "../ref", pattern=glob2rx("oeHBC_*exclude.Rda"))

reg3g_list_all = NULL;
omp_cyp_list_all = NULL;
csf1r_eGFP_list_all = NULL;

sapply(exfiles, function(exfile){
  lists <- load(file.path("../ref",exfile))
  reg3g_list_all <<- union(reg3g_list_all,reg3g_list)
  omp_cyp_list_all <<- union(omp_cyp_list_all, omp_cyp_list)
  csf1r_eGFP_list_all <<- union(csf1r_eGFP_list_all, csf1r_eGFP_list)
  rm(list=lists)
})
save(reg3g_list_all, omp_cyp_list_all, csf1r_eGFP_list_all, file="oeHBC_exclude.Rda")
