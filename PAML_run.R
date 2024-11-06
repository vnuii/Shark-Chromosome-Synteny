# run yn00
library(parallel)
library(dplyr)
ctl_names <- list.files("/home/vonui/odp_shark/ka_ks/PAML/all_InAndOut", pattern = ".ctl", full.names = T, recursive = T)
mclapply(ctl_names, function(i) {
  system(paste0("/home/vonui/anaconda3/envs/paml/bin/yn00 ", i))
}, mc.cores = 20)
