# script prepping occupancy model output data for processing

# Load libraries ----------------------------------------------------------

# check libraries carfully. this is especially important on the server/cluster
#.libPaths("/home/tietje/mamba/envs/r453/lib/R/library") 

rm(list=(ls())) # clear workspace
library(data.table)
library(terra)
library(sf)
library(exactextractr)

source("99_functions.R")

did = 5
if(did==6){scalID=2}else{scalID=1}



# Load & process data --------------------------------------------------------------------

# chose a did: 5(CZ), 6(NY), 17(NZ), 26(EU)

files = dir(paste0('data/occ_', as.character(did)), full.names = T)

tmp = list()
for(i in 1:length(files)){
  org = readRDS(files[i])
  tmp2 = list()
  for(l in 1:length(org)){
    scientificName = unlist(org[[l]][[2]]$sp_name)
    raw.occ = unlist(org[[l]][[2]]$raw.occ)
    tmp_int = cbind(data.table(scientificName, raw.occ), org[[l]][[1]])
    setDT(tmp_int)
    tmp2[[l]] = tmp_int
  }
  
  tmp[[i]] = rbindlist(tmp2) 
  tmp[[i]]$scientificName = gsub(pattern = " ", "_", tmp[[i]]$scientificName)
  tmp[[i]]$datasetID = did
  tmp[[i]]$scalingID = 1
  tmp[[i]]$endYear = i
}
rm(org)

dat = rbindlist(tmp)
rm(tmp, tmp2)
gc()





## Check for NAs 
apply(dat, 2, function(x){any(is.na(x))})
lunique(dat$scientificName)

## Remove NAs
dat = na.omit(dat)
apply(dat, 2, function(x){any(is.na(x))})
lunique(dat$scientificName)


# subset to species that occur in first and last time periods

## count years with data
times = range(dat$endYear)
dat[, alltime:=all(times %in% endYear), by=.(scientificName)]
tmp = unique(dat[,.(scientificName, alltime)])
table(tmp$alltime)

# not all species occur in first and last sampling period: remove
dat <- dat[alltime==TRUE, ]
lunique(dat$scientificName)

# more sanity checks
## any species that have only NAs? should be all FALSE
table(tapply(dat$mean.psi, dat$scientificName, function(x){all(is.na(x))}))
# and NA values in mean.psi values (occupancy values)? should be FALSE
any(is.na(dat$mean.psi))


## Subset to first and last period. releveant for quality testing
dat = dat[endYear %in% times,]



# -- Visual inspection --------------------------------------------

# load flagged species and exclude them
flagged_sp = gsub(".png", "", dir(paste0("maps/atlas=", did ,"_flagged_visually/")))
flag = data.table(species = gsub(".*[0-9]_", "", flagged_sp),
                  atlas = as.numeric(gsub(".*=|_[A-Z].*", "", flagged_sp)))
dat = dat[!scientificName %in% flag$species,]




# map quality check  -------------------------------------------------------

# check for high sd, mean.psi being off from raw.occ and lack of value variation
# (those are very stable species which make correlation unrealible)

check = unique(dat[, .(flos_flag=ifelse(sd.psi >= 0.20 & mean.psi > raw.occ * 10, "flag", "ok"),
                       psi.range=diff(range(mean.psi))),
                   .(scientificName, endYear, datasetID)])

sp <- check[flos_flag=="flag", scientificName]
sp2 <- check[psi.range<0.05, scientificName]
s = unique(c(sp,sp2))
dat <- dat[!scientificName %in% s, ]

lunique(dat$scientificName)



##### SAVE ###########

dat_s <- dat[, c('scientificName', 'siteID', 'pres.abs', 'raw.occ', 'mean.psi', 'sd.psi', 
        'geometry', 'datasetID', 'scalingID', 'endYear')]
saveRDS(dat_s, 
        paste0("data/processed_occupancy_for_CFs_2_",did,".rds"))














# Add elevation to grids --------------------------------------------------

for(did in c(5,6,17,26)){
  
  if(did==6){scalID=2}else{scalID=1}
  
  dat = readRDS(paste0("data/processed_occupancy_for_CFs_2_",did,".rds"))
  grid = st_read(paste0("data/all_scales_atlas_",did,".gpkg"))
  grid = grid[grid$scalingID==scalID,]
  grid = grid[grid$site %in% dat$siteID, ]
  
  ele <- terra::rast("data/counterfactuals/env/elevation/wc2.1_30s/wc2.1_30s_elev.tif")
  ele <- terra::crop(ele, grid)
  elevation <- exactextractr::exact_extract(ele, grid, fun="mean")
  grid = cbind(grid, elevation)
  write_sf(grid, paste0("data/one_scale_atlas_with_elevation",did,".gpkg"), )
  
  cat(did, "\r")
}





