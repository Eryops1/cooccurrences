# climate models


# Load libraries ----------------------------------------------------------

rm(list=(ls())) # clear workspace
library(data.table)
library(terra)
library(sf)
library(geodata)
library(ggplot2)
library(caret)
library(exactextractr)
library(parallel)
library(gbm)


source("99_functions.R")




# Load data --------------------------------------------------------------------


dataset_id = 26 # adjust as needed: 5,6,17,26

if(dataset_id==6){scalID=2}else{scalID=1}

dat = readRDS(paste0("data/processed_occupancy_for_CFs_2_",dataset_id,".rds"))
files = dir("data/environment/climate/annual", full.names = T, pattern = paste0("^" ,dataset_id, ".*(T1|T2)\\.tif"))
clim <- rast(files)
nam = gsub(paste0(".*annual/",dataset_id,"_|\\.tif"), "", files)
names(clim) = nam

# land use
if(dataset_id==5){years=c("1987|2016")}
if(dataset_id==6){years=c("1983|2003")}
if(dataset_id==17){years=c("1974|2002")}
if(dataset_id==26){years=c("1984|2015")}

files = dir("data/landuse", full.names = T, pattern = years)
lu <- rast(files)









# Extract and attach to spatial grid ----------------------------------------------


grid = vect(paste0("data/all_scales_atlas_",dataset_id,".gpkg"))
grid = grid[grid$scalingID==scalID,]
grid = grid[grid$siteID %in% dat$siteID,]
crs(clim, proj=TRUE)
crs(grid, proj=TRUE)
crs(lu, proj=TRUE)

#grid = extract(clim, grid, bind=TRUE, fun=mean)
#we gonna need something fast here since EU is in the mix
exclim = exactextractr::exact_extract(clim, st_as_sf(grid), fun="mean")
lu = crop(lu, grid)

target_classes <- c(11, 22, 33, 44, 55, 66, 77)
calc_props <- function(values, coverage_fractions) {
  # `values` is a data frame here because landcover has 2 layers
  out <- lapply(names(values), function(lyr) {
    v <- values[[lyr]]
    denom <- sum(coverage_fractions[!is.na(v)])  # or sum(coverage_fractions) if you want NA cells counted in the denominator
    p <- sapply(target_classes, function(cl) {
      sum(coverage_fractions[v == cl], na.rm = TRUE) / denom
    })
    setNames(p, paste0(lyr, "_frac_", target_classes))
})
as.data.frame(as.list(unlist(out)))
}
names(lu) = c("T1", "T2")
exlu <- exact_extract(lu, st_as_sf(grid), calc_props, force_df = TRUE, append_cols = "siteID")

names(exclim) = gsub("mean\\.", "", names(exclim))
grid = cbind(grid, exclim)
grid = cbind(grid, exlu)

# check out some maps
grid_sf = st_as_sf(grid)
ggplot(grid_sf, aes(fill=T2_frac_66))+
  geom_sf(col=NA)+
  theme_void()

ggplot(grid_sf, aes(fill=MAT_T2))+
  geom_sf(col=NA)+
  theme_void()

grid_dt = as.data.table(grid)


saveRDS(grid_dt, paste0("data/environment/climate/env_input_", dataset_id, ".rds"))








# Model species climate models --------------------------------------------


## GBM -----------
# I am trying a simple GBM first and truncate results. betareg is more
# technically correct, but values rarely exceed much above 1 and under 0. for
# our purpose, this is sufficient

dat = readRDS(paste0("data/processed_occupancy_for_CFs_2_",dataset_id,".rds"))
dat = merge(dat, grid_dt[, .(siteID, MAT_T1, MAT_T2, TAP_T1, TAP_T2, SEA_T_T1, 
                             SEA_T_T2, SEA_P_T1, SEA_P_T2,
                             T1_frac_11, T1_frac_22, T1_frac_33, T1_frac_44, T1_frac_55,
                             T1_frac_66, T1_frac_77, T2_frac_11, T2_frac_22, T2_frac_33,
                             T2_frac_44, T2_frac_55, T2_frac_66, T2_frac_77)], all.x=TRUE)


# standard parameters
fitControl_tune <- trainControl(
  method = "repeatedcv", 
  number = 10, 
  repeats = 10)  # full search, once per species

fitControl_fit  <- trainControl(
  method = "cv",
  number = 5) # light check, per simulation

ncores <- 20
species <- sort(unique(dat$scientificName))

for (i in seq_along(species)) {
  
  sim <- readRDS(paste0("data/uncertainty/simulated_", species[i], "_", dataset_id, ".rds"))
  
  dat_t1 <- dat[scientificName == species[i] & endYear == min(endYear)]
  setnames(dat_t1,
           old = c('MAT_T1','TAP_T1','SEA_T_T1','SEA_P_T1',
                   "T1_frac_11","T1_frac_22","T1_frac_33","T1_frac_44",
                   "T1_frac_55","T1_frac_66","T1_frac_77"),
           new = c('MAT','TAP','SEA_T','SEA_P','frac11','frac22','frac33','frac44',
                   'frac55','frac66','frac77'))
  
  dat_t2 <- dat[scientificName == species[i] & endYear == max(endYear)]
  setnames(dat_t2,
           old = c('MAT_T2','TAP_T2','SEA_T_T2','SEA_P_T2',
                   "T2_frac_11","T2_frac_22","T2_frac_33","T2_frac_44",
                   "T2_frac_55","T2_frac_66","T2_frac_77"),
           new = c('MAT','TAP','SEA_T','SEA_P','frac11','frac22','frac33','frac44',
                   'frac55','frac66','frac77'))
  
  ## ---- ONE-TIME TUNING PASS per species ----
  ## use the first simulated draw (or swap for the observed mean.psi if you have one)
  ## to find good hyperparameters via the full repeated CV search
  dat_t1_tune <- copy(dat_t1)
  dat_t1_tune$mean.psi <- sim[[1]][, 1]
  
  set.seed(825)
  tuneFit <- train(mean.psi ~ MAT+TAP+SEA_T+SEA_P+frac11+frac22+frac33+frac44+
                     frac55+frac66+frac77,
                   data = dat_t1_tune,
                   method = "gbm",
                   trControl = fitControl_tune,
                   verbose = FALSE)
  
  whichTwoPct <- tolerance(tuneFit$results, metric = "Rsquared", tol = 2, maximize = TRUE)
  bestGrid <- tuneFit$results[whichTwoPct, c("n.trees", "interaction.depth", "shrinkage", "n.minobsinnode")]
  
  ## ---- 100 simulations, hyperparameters fixed, light CV just for an R² ----
  
  # sim_mat = matrix(nrow=length(dat_t1_tune$mean.psi), ncol=100)
  # mod_out = list()
  # 
  sim_res <- mclapply(1:100, function(s) {
    
    tryCatch({
      tmp_t1 <- copy(dat_t1)
      tmp_t1$mean.psi <- sim[[1]][, s]
      
      tmp_t2 <- copy(dat_t2)
      tmp_t2$mean.psi <- last(sim)[, s]
      
      set.seed(825)
      gbmFit1 <- train(mean.psi ~ MAT+TAP+SEA_T+SEA_P+frac11+frac22+frac33+frac44+
                         frac55+frac66+frac77,
                       data = tmp_t1,
                       method = "gbm",
                       trControl = fitControl_fit,
                       tuneGrid = bestGrid,   # <- no grid search, single combo
                       verbose = FALSE)
      
      clim_occ <- predict(gbmFit1, newdata = tmp_t2)
      clim_occ <- pmax(0, pmin(1, clim_occ))
      #tmp_t2[, mean.psi_clim := clim_occ]
      #tmp_t2[, gbm_fit_Rsquared := gbmFit1$results$Rsquared]  # only one row now, no tolerance() needed
      #tmp_t2[, iteration := s]
      #sim_mat[, s] = clim_occ
      fit = caret::varImp(gbmFit1)#,
      gbm_fit_Rsquared = gbmFit1$results$Rsquared
           # res = tmp_t2[, .(siteID, scientificName, pres.abs, raw.occ, mean.psi, sd.psi,
           #                  geometry, datasetID, scalingID, endYear,
           #                  mean.psi_clim, gbm_fit_Rsquared, iteration)])
    }, error = function(e) e)
    list(clim_sim = clim_occ, fit = fit, gbm_fit_Rsquared=gbm_fit_Rsquared)
  }, mc.cores = ncores)
  
  mod_stats   <- lapply(sim_res, `[[`, "fit")
  mod_res <- lapply(sim_res, `[[`, "clim_sim")
  mod_R2 <- lapply(sim_res, `[[`, "gbm_fit_Rsquared")
  mod_res   <- matrix(unlist(mod_res), unique(lengths(mod_res)))
  
  saveRDS(list(mod_res, mod_stats, unlist(mod_R2)),
          paste0("data/counterfactuals2/c+l_sim_T2_species_", species[i], "_", dataset_id, ".rds"))
  
  cat("done:", species[i], "(", i, "/", length(species), ")\n")
}








