# Script: 01_get_co-occurrence_estimates.R
# Author: Melanie Tietje
# Email: tietje@fzp.czu.cz
# GitHub: @Eryops1
# Last Modified: 2025-10-27
# Purpose: Loads raw presence absence and occupancy probability data, saves the
#          processed data for analysis.
# Outputs: Four .rds objects named e.g. "more_lists_atlas=5_SES_cor_scales=1_2025-06-07.rds"
# Notes: groundhog will ensure the exact same R packages will be used and create
#        a library on first run, which might take a while.
#        Run this script for each atlas separately by adjusting the dataset_id variable

# Load libraries ----------------------------------------------------------

rm(list=(ls())) # clear workspace
library(groundhog)
pkgs <- c("data.table",
          "ggplot2",
          "terra",
          "sf",
          "parallel",
          "tictoc")
groundhog.library(pkgs, "2025-10-25")
theme_set(theme_classic())







# Load data --------------------------------------------------------------------

# chose a dataset_ids: 5(CZ), 6(NY), 17(NZ), 26(EU)
dataset_id = 26

files = dir(paste0('data/occ_', as.character(dataset_id)), full.names = T)

tmp = list()
for(i in 1:length(files)){
  # add scientificName
  org = readRDS(files[i])
  tmp2 = list()
  for(l in 1:length(org)){
    scientificName = unlist(org[[l]][[2]]$sp_name)
    tmp_int = cbind(org[[l]][[1]], scientificName)
    setDT(tmp_int)
    tmp2[[l]] = tmp_int
  }
  
  tmp[[i]] = rbindlist(tmp2) 
  tmp[[i]]$scientificName = gsub(pattern = " ", "_", tmp[[i]]$scientificName)
  tmp[[i]]$datasetID = dataset_id
  tmp[[i]]$scalingID = 1
  tmp[[i]]$endYear = i
}
rm(org)

dat = rbindlist(tmp)
rm(tmp)
gc()


# remove Russia and eastern Europe due to sampling effort differences (only relevant for EBBA)
if(dataset_id==26){
  grid = vect(paste0("data/all_scales_atlas_", dataset_id, ".gpkg"))
  grid = grid[grid$datasetID==26,]
  grid = grid[grid$scalingID==1,]

  all(dat$siteID %in% grid$siteID)
  plot(grid[grid$siteID %in% dat$siteID])
  # looks good?

  # grid = vect(paste0("data/all_scales_atlas_", dataset_id, ".gpkg"))
  # grid = grid[grid$datasetID==26,]
  # grid = grid[grid$scalingID==1,]
  # # build parts to remove
  # files = dir("data", pattern = "gadm41", full.names = T)
  # rem = lapply(files, vect)
  # rem = vect(rem)
  # # match projection
  # crs(rem, proj=T) == crs(grid, proj=T)
  # # intersect
  # tmp = terra::intersect(grid,rem)
  # tmp$remove = 1
  # grid_new = merge(grid, tmp, all.x=T)
  # grid = grid_new[is.na(grid_new$remove),]
}

# set scaling ID to 2 for the NY atlas which is the one processed for occupancy probabilities
if(dataset_id==6){
 dat$scalingID=2
}





# subset to species that occur in all time periods
# step 1: remove all NAs
dat = na.omit(dat)
# step2: count years with data
dat[, alltime:=length(unique(endYear)), by=.(scientificName)]
alltime_sp = unique(dat$scientificName[dat$alltime==length(files)])
dat <- dat[scientificName %in% alltime_sp,]

lunique = function(x){length(unique(x))}
end_years = sort(unique(dat$endYear))
scales = sort(unique(dat$scalingID))


# some sanity checks
## any species that have only NAs? should be all FALSE
table(tapply(dat$mean.psi, dat$scientificName, function(x){all(is.na(x))}))
# and NA values in mean.psi values (occupancy values)? should be FALSE
any(is.na(dat$mean.psi))


## plot map of a random species to check if everything looks ok
grid = vect(paste0("data/all_scales_atlas_", dataset_id, ".gpkg")) # load atlas grid
tmp = dat[datasetID==dataset_id & scientificName==unique(dat$scientificName)[2],] # select random species (in this case 2)
p1 = merge(grid, tmp, by=c("datasetID", "scalingID", "siteID"), all.x=TRUE) # merge into spatial grid

### presence absence data
ggplot(st_as_sf(p1[p1$scalingID==1,]), aes(fill=pres.abs))+
  geom_sf()+
  ggtitle(tmp$scientificName)
### occupancy probability data
ggplot(st_as_sf(p1[p1$scalingID==1,]), aes(fill=mean.psi))+
  geom_sf()+
  ggtitle(tmp$scientificName)












# Get co-occurrence measures ---------------------------------------------------

# we gonna get spearman correlation and C-scores (vegdist) 

# create species pairs 
nam = unique(dat$scientificName)
nam = nam[nam!=""] # no empties
pairs = combn(nam, 2, simplify = T) # matrix object with ncol= n pairs and row=2 (each column one pair)


# get grid of scale time combinations and chose species pairs
n = 100 # number of randomization for null distribution z score
s = unique(dat$scalingID)
t = seq(1:lunique(dat$endYear))
combs = expand.grid(s,t) # a dataframe with all combinations of scales and sampling periods
pairssub = pairs[,1:20] # subset of pairs for testing

## Define function for parallelization ####
get_cooccurrences <- function(dat, pairs, n, l) {
  library(data.table)
  # dat: data.table with siteID, pres.abs, effort, occupancy measures
  #       (mean.psi) and scientificName of species 
  # pairs: matrix object of species pairs, dimension: 2 rows, cols = number of 
  #         pairs  
  # n: iterations for null model (z-scores)
  # l: go through lines in combs object indicating scale and sampling period 
  # 
  scale = as.character(combs[l,1])
  time = combs[l,2]
  tmp = dat[endYear == time & scalingID == scale, ]
  res_int = list()
  
  for (i in 1:ncol(pairs)) { # for every species pair
    
    # get the data for species pair: sp1, sp2
    sp1 = sf::st_drop_geometry(tmp[scientificName %in% pairs[, i][1]])
    data.table::setnames(sp1, old=c('mean.psi', 'pres.abs'), new=c('mean.psi_sp1', 'pres.abs_sp1'))
    sp2 = sf::st_drop_geometry(tmp[scientificName %in% pairs[, i][2]])
    data.table::setnames(sp2, old=c('mean.psi', 'pres.abs'), new=c('mean.psi_sp2', 'pres.abs_sp2'))
    tmp2 = merge(sp1[, .(mean.psi_sp1, pres.abs_sp1, siteID)], sp2[, .(mean.psi_sp2, pres.abs_sp2, siteID)], all=T)
    
    # test for spatial overlap in presence grids
    ov = rowSums(tmp2[, .(pres.abs_sp1, pres.abs_sp2)])
    overlap = ifelse(all(ov<2), FALSE, TRUE)
    
    # get Spearman correlation and C-score
    cor_obs = cor(tmp2$mean.psi_sp1, tmp2$mean.psi_sp2, method='s')
    cor_obs_p = cor.test(tmp2$mean.psi_sp1, tmp2$mean.psi_sp2, method='s')$p.value
    c_obs = bipartite::C.score(cbind(tmp2$pres.abs_sp1, tmp2$pres.abs_sp2), normalise = TRUE)

    # null shuffle for Z-scores
    cor_null = rep(NA, n)
    cscore_null = rep(NA, n)
    for (j in 1:n) { # get n randomized co-occurrence measures
      rand = as.data.frame(tmp2) # turn data.frame because the sample function works different in DT objects
      rand = as.data.frame(apply(tmp2, 2, sample)) # shuffle all columns
      cor_null[j] = cor(rand$mean.psi_sp1, rand$mean.psi_sp2, method='s')
      cscore_null[j] = bipartite::C.score(cbind(rand$pres.abs_sp1, rand$pres.abs_sp2), normalise = TRUE)
    }
    
    # calculate Z-scores from the null distributions
    ses_cor = (cor_obs - mean(cor_null)) / sd(cor_null)
    ses_c = (c_obs - mean(cscore_null)) / sd(cscore_null)
    
    # collect results in data.table object
    res = data.table::data.table(
      species_pair = paste(pairs[, i], collapse = "|"),
      cor_obs = cor_obs,
      cor_obs_p = cor_obs_p,
      cor_ses = ses_cor,
      c_obs = c_obs,
      c_ses = ses_c,
      time_bin = paste0("T", time),
      scaleID = paste0("S", scale),
      dataset_id = dataset_id,
      overlap = overlap
    )
    
    # collect results data.table and the null distributions estimates in list object
    res_int[[i]] = list(res, cor_null, cscore_null)
    
    # print process. might not work in parallel as expected
    cat(paste0("time bin=", time, ", scale=", scale, ", species pair=", i), "\r")
  }
  return(res_int) # return results
}




## run ----

if(Sys.info()['sysname']=="Linux"){
  ncores = ifelse(nrow(combs)>16, 16, nrow(combs)) # cap cores at 16
  tic()
  res <- mclapply(1:nrow(combs), function(l) {
    get_cooccurrences(dat = dat, pairs = pairs, n = n, l)
  }, mc.cores = ncores)
  toc() # ca 17 minutes for CZ scale 1
}

if(Sys.info()['sysname']=="Windows"){
nclust = ifelse(nrow(combs)>16, 16, nrow(combs)) # cap cores at 16
cl <- parallel::makeCluster(nclust)
parallel::clusterExport(cl, 
                        varlist=list("dat", 
                                     "combs", 
                                     "pairs", 
                                     "pairssub", 
                                     "get_cooccurrences", 
                                     "dataset_id"))
tic()
res <- parLapply(cl, 1:nrow(combs), function(l) {
  get_cooccurrences(dat=dat, pairs = pairs, n = 100, l)
})
toc()
stopCluster(cl)
}




## save ----
outfile = paste0("data/more_lists_atlas=", dataset_id, "_SES_cor_scales=", paste0(s, collapse = "_"), "_", Sys.Date(), ".rds")
saveRDS(res, outfile)


