# Script: 02_get_changes_in_co-occurrences.R
# Author: Melanie Tietje
# Email: tietje@fzp.czu.cz
# GitHub: @Eryops1
# Last Modified: 2025-11-05
# Purpose: Loads co-occurrence estimates, gets changes between sampling periods
#          incl null model estimates, saves the processed data for analysis.
# Output: .rds object named "data/processed_spass.rds"
# Notes: groundhog will ensure the exact same R packages will be used and create
#        a library on first run, which might take a while
#        Runs parts in parallel, tested for linux should work 

# Load libraries ----------------------------------------------------------

rm(list=(ls())) # clear workspace
library(groundhog)
pkgs <- c("data.table",
          "parallel",
          "tictoc")
groundhog.library(pkgs, "2025-10-25")



# functions --------------------------------------------------------------
sort_pairs_vectorized <- function(x) {
  # a function that sorts species pair names string into alphabetic order
  sapply(x, function(item) {
    spl <- unlist(strsplit(item, split = "|", fixed = TRUE))
    paste0(sort(spl), collapse = "|")
  })
}

process_change <- function(a, org) {
  res_int <- vector("list", length(org[[a]][[1]]))
  
  for (s in seq_along(org[[a]][[1]])) {
    nulldat <- data.table()
    
    for (i in seq_along(org[[a]])) {
      nulldat <- rbind(nulldat, cbind(
        org[[a]][[i]][[s]][[1]],
        org[[a]][[i]][[s]][[2]],
        org[[a]][[i]][[s]][[3]]
      ))
    }
    
    setnames(nulldat, c("V2", "V3"), c("cor_null", "c_null"))
    nulldat$key <- rep(1:100, length(org[[a]]))
    
    nulldat[, delta_cor_null := cor_null[time_bin == tail(time_bin, 1)] - cor_null[time_bin == "T1"], by = .(key, scaleID)]
    nulldat[, delta_cor_obs := cor_obs[time_bin == tail(time_bin, 1)] - cor_obs[time_bin == "T1"], by = .(key, scaleID)]
    nulldat[, delta_cor_obs_z := (delta_cor_obs - mean(delta_cor_null)) / sd(delta_cor_null), by = .(scaleID)]
    nulldat[, delta_cor_ses := cor_ses[time_bin == tail(time_bin, 1)] - cor_ses[time_bin == "T1"], by = .(key, scaleID)]
    
    nulldat[, delta_c_null := c_null[time_bin == tail(time_bin, 1)] - c_null[time_bin == "T1"], by = .(key, scaleID)]
    nulldat[, delta_c_obs := c_obs[time_bin == tail(time_bin, 1)] - c_obs[time_bin == "T1"], by = .(key, scaleID)]
    nulldat[, delta_c_obs_z := (delta_c_obs - mean(delta_c_null)) / sd(delta_c_null), by = .(scaleID)]
    nulldat[, delta_c_ses := c_ses[time_bin == tail(time_bin, 1)] - c_ses[time_bin == "T1"], by = .(key, scaleID)]
    
    res_int[[s]] <- unique(nulldat[, .(
      species_pair, scaleID, time_bin,
      cor_obs, cor_obs_p, cor_ses, delta_cor_obs,
      delta_cor_obs_z, delta_cor_ses,
      c_obs, c_ses, delta_c_obs,
      delta_c_obs_z, delta_c_ses, overlap
    )])
  }
  
  return(rbindlist(res_int))
}













# Read data ---------------------------------------------------------------


files1 = dir("data", pattern = "cor_", full.names = TRUE)
files2 = files1[grep("chunk", files1)]

# little list magic code for the atlas 26 chunks
org1 = lapply(files2, readRDS)
o1 = c(org1[[1]][], org1[[2]][])
o2 = c(org1[[3]][], org1[[4]][])
org1 = list(o1,o2)

org2 = lapply(files1[!grepl("chunk", files1)], readRDS)
# combine
org = c(list(org1), org2)
rm(org1, org2)


# attach atlas name (double check)
atlas = c("26", "17", "5", "6")
names(org) =atlas



# Process data ---------------------------------------------------------------

# should take about 13 minutes on 4 cores

if(Sys.info()['sysname']=="Linux"){
  ncores = ifelse(length(atlas)>16, 16, length(atlas)) # cap cores at 16, just in case
  tic()
  res <- mclapply(1:length(atlas), function(a) {
    process_change(a, org)
  }, mc.cores = ncores)
  toc() 
}

if(Sys.info()['sysname']=="Windows"){
  nclust = ifelse(length(atlas)>16, 16, length(atlas)) # cap cores at 16
  cl <- parallel::makeCluster(nclust)
  parallel::clusterExport(cl, 
                          varlist=list("org", 
                                       "atlas", 
                                       "process_change"))
  tic()
  res <- parLapply(cl, 1:length(atlas), function(a) {
    process_change(a, org)
  })
  toc()
  stopCluster(cl)
}

names(res) = names(org)
res = rbindlist(res, idcol = "dataset_id")
res$dataset_id = as.numeric(res$dataset_id)

# number of species pairs per atlas
tapply(res$species_pair, res$dataset_id, function(x){length(unique(x))})


# SAVE ----
saveRDS(res, "data/processed_spass.rds")







