# checking status of simulations. test maps, counts, etc
library(data.table)
library(ggplot2)
library(cowplot)









spa <- data.table(total = rep(NA,4),
                  atlas = c("5", "6", "17", "26"))
nam = list()
for(i in 1:length(spa$atlas)){
  tmp <- readRDS(paste0("data/processed_occupancy_for_CFs_2_",spa$atlas[i],".rds"))
  spa$total[i] <- length(unique(tmp$scientificName))
  nam[[i]] = data.table(scientificName = sort(unique(tmp$scientificName)))
}
spa
names(nam) = spa$atlas
nam = rbindlist(nam, use.names = T, idcol = "dataset_id")
nam[, nam_id:=paste(scientificName, dataset_id, sep="_"),]




# how many are done -------------------------------------------------------

local_dir <- "data/counterfactuals2/" 

## CF_IM
files = dir(local_dir, pattern="simulated_[A-Z].*", full.names = T)
files = files[grep(paste0(nam$nam_id, collapse = "|"), files)] # select for target species
tmp = regmatches(files, gregexpr("simulated_[A-Z].*_[0-9]{1,2}", files))
tmp = gsub("simulated_[A-Z].*_[a-z].*_", "", tmp)
tmp <- as.data.table(table(tmp), keep.rownames = T)
setnames(tmp, old=c("tmp", "N"), new=c("atlas", "CF_IM"))
fin <- merge(spa, tmp, by="atlas")

## CF_IM_elev
files = dir(local_dir, pattern="simulated_elevation_[A-Z].*")
files = files[grep(paste0(nam$nam_id, collapse = "|"), files)] # select for target species
tmp = regmatches(files, gregexpr("simulated_elevation_[A-Z].*_[0-9]{1,2}", files))
tmp = gsub("simulated_elevation_[A-Z].*_[a-z].*_", "", tmp)
tmp <- as.data.table(table(tmp), keep.rownames = T)
setnames(tmp, old=c("tmp", "N"), new=c("atlas", "CF_IM_elev"))
fin <- merge(fin, tmp, by="atlas", all.x=T)

## CF_IM_env
files = dir(local_dir, pattern="simulated_env_[A-Z].*")
files = files[grep(paste0(nam$nam_id, collapse = "|"), files)] # select for target species
tmp = regmatches(files, gregexpr("simulated_env_[A-Z].*_[0-9]{1,2}", files))
tmp = gsub("simulated_env_[A-Z].*_[a-z].*_", "", tmp)
tmp <- as.data.table(table(tmp), keep.rownames = T)
setnames(tmp, old=c("tmp", "N"), new=c("atlas", "CF_IM_env"))
fin <- merge(fin, tmp, by="atlas", all.x=T)

## CF_climate (this is clim + lu but i keep calling it that for simpl)
files = dir(local_dir, pattern="^c\\+l*")
files = files[grep(paste0(nam$nam_id, collapse = "|"), files)] # select for target species
tmp = unlist(regmatches(files, gregexpr("[0-9]{1,2}", files)))
tmp <- as.data.table(table(tmp), keep.rownames = T)
setnames(tmp, old=c("tmp", "N"), new=c("atlas", "C+L"))
fin <- merge(fin, tmp, by="atlas", all.x=T)


fin












# check sim dimensions -------------------------------------------------------

## CF_IM
files = dir(local_dir, pattern="simulated_[A-Z].*", full.names = T)
files = files[grep(paste0(nam$nam_id, collapse = "|"), files)] # select for target species

dims = list()
for(i in 1:length(spa$atlas)){
  tmp = lapply(files[grep(paste0("_", spa$atlas[i]), files)], readRDS)
  dims[[i]] = unique(matrix(unlist(lapply(tmp, dim)), ncol=2, byrow = T))
}
dims


## CF_IM_elev
files = dir(local_dir, pattern="simulated_elevation_[A-Z].*", full.names = T)
files = files[grep(paste0(nam$nam_id, collapse = "|"), files)] # select for target species

dims = list()
for(i in 1:length(spa$atlas)){
  tmp = lapply(files[grep(paste0("_", spa$atlas[i]), files)], readRDS)
  dims[[i]] = unique(matrix(unlist(lapply(tmp, dim)), ncol=2, byrow = T))
}
dims

## CF_IM_env
files = dir(local_dir, pattern="simulated_env_[A-Z].*", full.names = T)
files = files[grep(paste0(nam$nam_id, collapse = "|"), files)] # select for target species

dims = list()
for(i in 1:length(spa$atlas)){
  tmp = lapply(files[grep(paste0("_", spa$atlas[i]), files)], readRDS)
  dims[[i]] = unique(matrix(unlist(lapply(tmp, dim)), ncol=2, byrow = T))
}
dims


# all simulations within atlases are of the same and correct dimensions









# check number of succesful sims ------------------------------------------

# how many NA cols?

cfs = c("CF_IM","CF_IM_elev","CF_IM_env")
pats = c("simulated_[A-Z].*",
         "simulated_elevation_[A-Z].*", 
         "simulated_env_[A-Z].*")

nsf = list()
for(cf in 1:length(cfs)){
  p = pats[cf]
  files = dir(local_dir, pattern=p, full.names = T)
  files = files[grep(paste0(nam$nam_id, collapse = "|"), files)] # select for target species
  int_list = list()
  
  for(i in 1:length(spa$atlas)){
    filessub = files[grep(paste0("_", spa$atlas[i]), files)]
    tmp = lapply(filessub, readRDS)
    tmp2 = lapply(tmp,function(x){
      apply(x, 2, function(x)any(is.na(x)))
    })
    int_list[[i]] = data.table(atlas = spa$atlas[i],
                               treatment = cfs[cf],
                          scientificName = regmatches(filessub, gregexpr("[A-Z][a-z].*\\_[a-z].*\\_", filessub)),
                          n_complete_sims = unlist(lapply(tmp2, function(x)length(which(x!=T)))))
  }
  
  nsf[[cf]] = rbindlist(int_list)
  cat(cf, "\r")
}

nsf = rbindlist(nsf)
nsf[ order(n_complete_sims, decreasing = F), ]



saveRDS(nsf, "data/n_complete_sims.rds")










# maps --------------------------------------------------------------------

# make set of maps for example species from the original and simulations

# select dataset
dataset_id = 5
if(dataset_id==6){scalID=2}else{scalID=1}

dat = readRDS(paste0("data/processed_occupancy_for_CFs_2_",dataset_id,".rds"))
grid = vect(paste0("data/all_scales_atlas_",dataset_id,".gpkg"))
grid = grid[grid$scalingID==scalID,]


# select species
species = sort(unique(dat$scientificName))
sp = species[5]
sub = dat[scientificName ==sp,]




# read simulations
files <- dir("data/counterfactuals2", paste0(sp, "_", dataset_id, ".rds"), full.names = T)
sim = lapply(files, readRDS)
names(sim) = gsub("data/counterfactuals2/|_[A-Z].*", "", files)

# read uncertainty
fules <- dir("data/uncertainty", paste0("simulated_", sp, "_", dataset_id, ".rds"), full.names = T)
un = readRDS(fules)


# grid with uncertainty
grid = grid[grid$siteID %in% sub$siteID,]
un1 = as.data.table(un[[1]])
names(un1) = paste0("uncertainty_T1_", names(un1))
un2 = as.data.table(last(un))
names(un2) = paste0("uncertainty_T2_", names(un2))
grid_uncertainty = cbind(grid, un1, un2)

# grid with average model values
grid_mean_model = merge(grid, sub, all=T)

# grid with simulations
cf_im = as.data.table(sim$simulated)
names(cf_im) = paste0("CF_IM_", names(cf_im))
cf_im_elev = as.data.table(sim$simulated_elevation)
names(cf_im_elev) = paste0("CF_IM_elev_", names(cf_im_elev))
cf_im_env = as.data.table(sim$simulated_env)
names(cf_im_env) = paste0("CF_IM_env_", names(cf_im_env))
cf_im_clim = as.data.table(sim$`c+l_sim`[[1]])
names(cf_im_clim) = paste0("CF_clim_", names(cf_im_clim))

grid_sims = cbind(grid, cf_im, cf_im_elev, cf_im_env, cf_im_clim)



p_mean_model = st_as_sf(grid_mean_model)
p_uncertainty = st_as_sf(grid_uncertainty)
p_sims = st_as_sf(grid_sims)

plot_grid(nrow=3, rel_heights = c(1,1,0.7),
  # original maps
  plot_grid(ncol=2, #labels = "AUTO",
            ggplot(p_mean_model[p_mean_model$endYear==1,], aes(fill=mean.psi))+
              geom_sf(col=NA, show.legend = F)+
              theme_void()+
              ggtitle("Mean occupancy model, T1"),
            ggplot(p_mean_model[p_mean_model$endYear==max(p_mean_model$endYear),], aes(fill=mean.psi))+
              geom_sf(col=NA)+
              theme_void()+
              ggtitle("Mean occupancy model, T2")
            ),
  # uncertainty examples
  plot_grid(ncol=2,
            ggplot(p_uncertainty, aes(fill=uncertainty_T1_V1))+
              geom_sf(col=NA, show.legend = F)+
              theme_void()+
              ggtitle("Uncertainty map example, T1"),
            ggplot(p_uncertainty, aes(fill=uncertainty_T2_V1))+
              geom_sf(col=NA)+
              theme_void()+
              scale_fill_continuous("mean.psi")+
              ggtitle("Uncertainty map example, T2")
    ),
  # simulation examples (this is all T2 only)
  plot_grid(ncol=4,
    ggplot(p_sims, aes(fill=CF_clim_V1))+
      geom_sf(col=NA, show.legend = F)+
      theme_void()+
      ggtitle("Clim & LU prediction T2"),
    ggplot(p_sims, aes(fill=CF_IM_V1))+
      geom_sf(col=NA, show.legend = F)+
      theme_void()+
      ggtitle("CF_IM T2"),
    ggplot(p_sims, aes(fill=CF_IM_elev_V1))+
      geom_sf(col=NA, show.legend = F)+
      theme_void()+
      ggtitle("CF_IM_elevation T2"),
    ggplot(p_sims, aes(fill=CF_IM_env_V1))+
      geom_sf(col=NA, show.legend = F)+
      theme_void()+
      ggtitle("CF_IM_environment T2")
  )
)






