# counterfactuals to run on server on a lot of cores...

# Load libraries -----------------------------------------

.libPaths("/home/tietje/mamba/envs/r453/lib/R/library") 

rm(list=(ls())) # clear workspace
library(data.table)
library(terra)
library(sf)
library(parallel)
library(MASS)
library(fields)
library(sp)
library(RhpcBLASctl)
library(gstat)
#library(exactextractr) # pything lib missing, use terra nontheless


# restrict this to avoid BLAS thread oversubscription
blas_set_num_threads(4)
omp_set_num_threads(4)


print(paste("System time: ", Sys.time()))


# Bash input ###############################

(args <- commandArgs(trailingOnly=TRUE))

dataset_id    =as.numeric(args[1])
which_CF      =args[2]

nnull         =as.numeric(args[3])
keep_tol      =as.numeric(args[4])
target_n      =as.numeric(args[5])
batch_size    =as.numeric(args[6])
max_batches   =as.numeric(args[7])
keep_tol_elev =as.numeric(args[9])

ncores        =as.numeric(args[8])


# dataset_id=6
# nnull=500
# keep_tol=0.05
# target_n=1
# batch_size=500
# max_batches=10
# keep_tol_elev=0.1
# ncores=5



#################################################

if(dataset_id==6){scalID=2}else{scalID=1}


if(which_CF=="uncertainty"){
  # uncertainty ---------------------------------------------------
  # creates two lists for each species for resampled occupancy probablities
  
  dat = readRDS(paste0("data/processed_occupancy_for_CFs_2_",dataset_id,".rds"))
  species <- sort(unique(dat$scientificName))
  
  # uncertainty is different from counterfactuals later since we will be adding
  # uncertainty to both sampling periods. Therefore, we will estimate the
  # variogram on the go instead of reading the file
  
  grid = st_read(paste0("data/all_scales_atlas_",dataset_id,".gpkg"))
  grid = grid[grid$scalingID==scalID,]
  
  cent <- sf::st_centroid(unique(dat$geometry))
  xy <- sf::st_coordinates(cent)
  D <- as.matrix(dist(xy))
  
  # takes about 1h on 10 cores for CZ
  mclapply(seq_along(species), function(i){
    
    tmp <- dat[scientificName == species[i], ] # we are sampling both years!
    years = range(dat$endYear)
    
    tmp_year_res = list()
    
    for(y in years){
      
      tmp_year = tmp[endYear==y, ]
      # Estimate distance decay with variogram
      df <- data.table(x = xy[,1],y = xy[,2], occ = tmp_year$mean.psi)
      
      sp::coordinates(df) <- ~x+y # turn columns into actual coordinates
      vg <- gstat::variogram(occ ~ 1, df)
      fit <- gstat::fit.variogram(vg, model = vgm(model = "Exp"))
      range_obs <- fit$range
      
      # get covariance matrix
      R <- exp(-D / range_obs)
      Sigma <- outer(tmp_year$sd.psi, tmp_year$sd.psi) * R
      diag(Sigma) <- diag(Sigma) + 1e-6 # add tiny variation for numerical robustness
      
      occ_keep = NULL
      n_keep = 0
      batch = 1
      
      while(n_keep < target_n & batch <= max_batches){
        
        # multivariate sampling
        tmp_year_sim <- t(MASS::mvrnorm(n = batch_size, mu = tmp_year$mean.psi, Sigma = Sigma))
        
        ## variograms to test if range_param is within accepted tolerance
        tmpvg <- data.table(x = xy[,1], y = xy[,2], tmp_year_sim)
        sp::coordinates(tmpvg) <- ~x+y # turn into coordinates
        
        vg_sims <- data.table(
          iteration  = seq_len(ncol(tmp_year_sim)),
          warning    = FALSE,
          range_param = NA_real_
        )
        
        for(vgi in seq_len(ncol(tmp_year_sim))){
          subtmp <- tmpvg[,vgi]
          names(subtmp) <- "response"
          vg <- gstat::variogram(response ~ 1, subtmp)
          fit <- withCallingHandlers(
            gstat::fit.variogram(vg, model = gstat::vgm(model = "Exp")),
            warning = function(w){
              vg_sims$warning[vgi] <<- TRUE
              invokeRestart("muffleWarning")
            }
          )
          vg_sims$range_param[vgi] <- fit$range[[1]]
          
        }
        # check which ones to keep
        vg_sims[, keep:=abs(range_param - range_obs) / range_obs < keep_tol]
        keep_idx = which(vg_sims$keep)
        
        # check if long enough already
        if(length(keep_idx) > 0){
          occ_keep <- cbind(occ_keep, tmp_year_sim[, keep_idx, drop = FALSE])
          n_keep <- ncol(occ_keep)
        }
        batch <- batch + 1
      }
      
      # save collected simulations and move on to next year
      tmp_year_res[[y]] <- occ_keep[,1:target_n, drop=FALSE]
      
    }
    
    saveRDS(tmp_year_res, file = paste0("data/uncertainty/simulated_", species[i], "_", dataset_id, ".rds"))
  }, mc.cores = ncores)
}



if(which_CF=="CF_IM"){
  # CF IM ---------------------------------------------------
  
  dat = readRDS(paste0("data/processed_occupancy_for_CFs_2_",dataset_id,".rds"))
  
  species <- sort(unique(dat$scientificName))
  
  # check which are done, modify species accordingly
  sp_done = dir("data/counterfactuals2/", pattern = paste0("simulated_[A-Z].*_", dataset_id, "\\.rds"))
  sp_done = gsub(paste0("simulated_|_", dataset_id, "\\.rds"), "", sp_done)
  species = species[!species %in% sp_done]
  print(species)
  
  grid = st_read(paste0("data/all_scales_atlas_",dataset_id,".gpkg"))
  grid = grid[grid$scalingID==scalID,]
  
  cent <- sf::st_centroid(unique(dat$geometry))
  xy <- sf::st_coordinates(cent)
  D <- as.matrix(dist(xy))
  
  mclapply(seq_along(species), function(i){
    
    # read sim
    dat = readRDS(paste0("data/uncertainty/simulated_", species[i], "_", dataset_id, ".rds"))
    tmp = last(dat)
    cf = matrix(nrow=nrow(tmp), ncol = ncol(tmp), NA)
    
    for(s in 1:ncol(tmp)){
      # get variogram for each simulation (=column)
      df <- data.table(x = xy[,1],y = xy[,2], occ = tmp[,s])
      sp::coordinates(df) <- ~x+y # turn columns into actual coordinates
      vg <- gstat::variogram(occ ~ 1, df)
      vg_fit = gstat::fit.variogram(vg, model = vgm(model = "Exp"))
      range_obs = vg_fit$range[1]
      
      ## covariance matrix
      Sigma <- exp(-D / range_obs)
      diag(Sigma) <- diag(Sigma) + 1e-6
      
      # adjust with elevation
      
      occ_keep = NULL
      n_keep = 0
      batch = 1
      
      while(n_keep < target_n & batch <= max_batches){
        
        ## simulate one batch
        latent <- MASS::mvrnorm(n = batch_size, mu = rep(0, nrow(tmp)), Sigma = Sigma)
        
        ## rank reassignment
        occ_sorted = sort(tmp[,s])
        latent_ranks = apply(latent, 1, rank)
        occ_cf = apply(latent_ranks, 2, function(x) occ_sorted[x])
        
        ## variograms
        tmpvg <- data.table(x = xy[,1], y = xy[,2], occ_cf)
        sp::coordinates(tmpvg) <- ~x+y # turn into coordinates
        
        vg_sims <- data.table(
          iteration  = seq_len(ncol(occ_cf)),
          warning    = FALSE,
          range_param = NA_real_
        )
        
        pass_idx <- NA_integer_  # index of first column that passes, if any
        
        for(vgi in seq_len(ncol(occ_cf))){
          subtmp <- tmpvg[,vgi]
          names(subtmp) <- "response"
          vg <- gstat::variogram(response ~ 1, subtmp)
          fit <- withCallingHandlers(
            gstat::fit.variogram(vg, model = gstat::vgm(model = "Exp")),
            warning = function(w){
              vg_sims$warning[vgi] <<- TRUE
              invokeRestart("muffleWarning")
            }
          )
          vg_sims$range_param[vgi] <- fit$range[[1]]
          
          # stop as soon as we find the first passing column
          if(abs(vg_sims$range_param[vgi] - range_obs) / range_obs < keep_tol){
            pass_idx <- vgi
            break
          }
        }
        
        # check if long enough already
        if(!is.na(pass_idx)){
          occ_keep <- cbind(occ_keep, occ_cf[, pass_idx, drop = FALSE]) # save the column that actually passed
          n_keep <- ncol(occ_keep)
        }
        batch <- batch + 1
      }
      
      ## store the 1 simulation
      cf[, s] = occ_keep[,1:target_n, drop=FALSE]
      
    }
    
    saveRDS(cf, file = paste0("data/counterfactuals2/simulated_", species[i], "_", dataset_id, ".rds"))
    
  }, mc.cores = ncores)
}



if(which_CF=="CF_IM_elev"){
  # CF ELEV ---------------------------------------------
  
  dat = readRDS(paste0("data/processed_occupancy_for_CFs_2_",dataset_id,".rds"))
  
  species <- sort(unique(dat$scientificName))
  
  # check which are done, modify species accordingly
  sp_done = dir("data/counterfactuals2/", pattern = paste0("simulated_elevation_[A-Z].*_", dataset_id, "\\.rds"))
  sp_done = gsub(paste0("simulated_elevation_|_", dataset_id, "\\.rds"), "", sp_done)
  species = species[!species %in% sp_done]
  print(species)
  
  grid = st_read(paste0("data/one_scale_atlas_with_elevation",dataset_id,".gpkg"))
  grid = grid[grid$scalingID==scalID,]
  
  cent <- sf::st_centroid(unique(dat$geometry))
  xy <- sf::st_coordinates(cent)
  D <- as.matrix(dist(xy))
  
  # how influential should elevation be?
  a=0.95
  
  # scale elevation
  elev_z <- scale(grid$elevation)
  
  mclapply(seq_along(species), function(i){
    
    # read sim
    dat = readRDS(paste0("data/uncertainty/simulated_", species[i], "_", dataset_id, ".rds"))
    tmp = last(dat)
    cf = matrix(nrow=nrow(tmp), ncol = ncol(tmp), NA)
    
    for(s in 1:ncol(tmp)){
      # correlation species with elevation
      rho_elevation <- cor(tmp[, s], elev_z, method = "s")[1]
      
      # get variogram for each simulation (=column)
      df <- data.table(x = xy[,1],y = xy[,2], occ = tmp[,s])
      sp::coordinates(df) <- ~x+y # turn columns into actual coordinates
      vg <- gstat::variogram(occ ~ 1, df)
      vg_fit = gstat::fit.variogram(vg, model = vgm(model = "Exp"))
      range_obs = vg_fit$range[1]
      
      ## covariance matrix
      Sigma <- exp(-D / range_obs)
      diag(Sigma) <- diag(Sigma) + 1e-6
      
      occ_keep = NULL
      n_keep = 0
      batch = 1
      
      while(n_keep < target_n & batch <= max_batches){
        
        ## simulate one batch
        latent <- MASS::mvrnorm(n = batch_size, mu = rep(0, nrow(tmp)), Sigma = Sigma)
        
        ## ranking guassian field values, sorting occupancy values
        occ_sorted = sort(tmp[,s])
        latent_ranks = apply(latent, 1, rank)
        
        ## modify ranks with elevation
        alpha <- a * abs(rho_elevation)
        elev_rank <- rank(sign(rho_elevation) * elev_z)
        # NOTE: sign ensures ordering to be aligned with direction. low
        # ranks = low elevation, even with negative correaltion. this does not
        # matter of rho=positive, but for rh=negative
        score <- t(latent_ranks) +
          alpha * matrix(
            elev_rank,
            nrow=nnull,
            ncol=length(elev_rank),
            byrow=TRUE
          )
        score_ranks <- apply(score, 1, rank)
        occ_cf = apply(score_ranks, 2, function(x) occ_sorted[x])
        
        # test correlation
        occ_cf_elev_cor = apply(occ_cf, 2, function(x){cor(x, elev_z, method = "s")})
        
        ## variograms
        tmpvg <- data.table(x = xy[,1], y = xy[,2], occ_cf)
        sp::coordinates(tmpvg) <- ~x+y # turn into coordinates
        
        vg_sims <- data.table(
          iteration  = seq_len(ncol(occ_cf)),
          warning    = FALSE,
          range_param = NA_real_,
          cor_elevation = occ_cf_elev_cor
        )
        
        # check elevation correlation criterion FIRST (cheap) so we can skip
        # the expensive gstat fit for candidates that would be rejected anyway
        vg_sims[, keep_elev:=abs(cor_elevation - rho_elevation) < keep_tol_elev]
        vg_sims[, keep_sac:=FALSE]  # default; only set TRUE where we actually fit
        
        still_needed <- target_n - n_keep  # how many more acceptable columns needed
        
        for(vgi in seq_len(ncol(occ_cf))){
          
          # skip the variogram fit entirely if elevation criterion already fails
          if(!vg_sims$keep_elev[vgi]) next
          
          subtmp <- tmpvg[,vgi]
          names(subtmp) <- "response"
          vg <- gstat::variogram(response ~ 1, subtmp)
          fit <- withCallingHandlers(
            gstat::fit.variogram(vg, model = gstat::vgm(model = "Exp")),
            warning = function(w){
              vg_sims$warning[vgi] <<- TRUE
              invokeRestart("muffleWarning")
            }
          )
          vg_sims$range_param[vgi] <- fit$range[[1]]
          vg_sims$keep_sac[vgi] <- abs(vg_sims$range_param[vgi] - range_obs) / range_obs < keep_tol
          
          if(isTRUE(vg_sims$keep_sac[vgi])){
            still_needed <- still_needed - 1
            if(still_needed <= 0) break  # got enough from this batch, stop fitting
          }
        }
        
        # check which ones to keep, for range param and correlation with elevation
        keep_idx = which(vg_sims$keep_sac & vg_sims$keep_elev)
        
        # check if enough simulations meet requirements
        # check indexing!
        if(length(keep_idx) > 0){
          occ_keep <- cbind(occ_keep, occ_cf[, keep_idx, drop = FALSE])
          n_keep <- ncol(occ_keep)
        }
        batch <- batch + 1
      }
      
      ## store the 1 simulation (or NA if nothing can be found!)
      if(is.null(occ_keep)){
        cf[, s] = rep(NA, nrow(cf))
      }else{
        cf[, s] = occ_keep[,1:target_n, drop=FALSE] 
      }
      
      
    }
    saveRDS(cf, file = paste0("data/counterfactuals2/simulated_elevation_", species[i], "_", dataset_id, ".rds"))
    
  }, mc.cores = ncores)
}




if(which_CF=="CF_IM_env"){
  # CF ELEV+CLIM+LANDUSE = environment -------------------------------------
  
  dat = readRDS(paste0("data/processed_occupancy_for_CFs_2_",dataset_id,".rds"))
  
  grid = vect(paste0("data/one_scale_atlas_with_elevation",dataset_id,".gpkg"))
  grid = grid[grid$scalingID==scalID,]
  
  species <- sort(unique(dat$scientificName))
  
  # check which are done, modify species accordingly
  sp_done = dir("data/counterfactuals2/", pattern = paste0("^simulated_env_.*_", dataset_id, "\\.rds"))
  sp_done = gsub(paste0("simulated_env_|_", dataset_id, "\\.rds"), "", sp_done)
  species = species[!species %in% sp_done]
  print(species)
  
  cent <- sf::st_centroid(unique(dat$geometry))
  xy <- sf::st_coordinates(cent)
  D <- as.matrix(dist(xy))
  
  # how influential should elevation be?
  a=0.95
  # scale elevation
  elev_z <- scale(grid$elevation)
  
  # how influential should environment be?
  a_clim=0.95
  
  mclapply(seq_along(species), function(i){

    # read uncertainty simulation
    dat = readRDS(paste0("data/uncertainty/simulated_", species[i], "_", dataset_id, ".rds"))
    tmp = last(dat) # sampling the last year

    # read environmental model simulation
    env_mod = readRDS(paste0("data/counterfactuals2/c+l_sim_T2_species_", species[i], "_", dataset_id, ".rds"))
    env = env_mod[[1]] # environmental model matrix, x by 100 columns. we use that as a guiding map

    cf = matrix(nrow=nrow(tmp), ncol = ncol(tmp), NA) # same dims as the input
    
    for(s in 1:ncol(tmp)){
    # correlation species with elevation
      rho_elevation <- cor(tmp[, s], elev_z, method = "s")[1]
      
      # get variogram for this simulation
      df <- data.table(x = xy[,1],y = xy[,2], occ = tmp[,s])
      sp::coordinates(df) <- ~x+y # turn columns into actual coordinates
      vg <- gstat::variogram(occ ~ 1, df)
      vg_fit = gstat::fit.variogram(vg, model = vgm(model = "Exp"))
      range_obs = vg_fit$range[1]
       
     ## covariance matrix
     Sigma <- exp(-D / range_obs)
     diag(Sigma) <- diag(Sigma) + 1e-6
    
     occ_keep = NULL
     n_keep = 0
     batch = 1
    
    while(n_keep < target_n & batch <= max_batches){
      
      ## simulate one batch
      latent <- MASS::mvrnorm(n = batch_size, mu = rep(0, nrow(tmp)), Sigma = Sigma)
      
      ## ranking guassian field values, sorting occupancy values
      occ_sorted = sort(tmp[,s])
      latent_ranks = apply(latent, 1, rank)
      
      ## modify ranks with elevation + climate + landuse, according to R2! 
      alpha <- a * abs(rho_elevation)
      alpha_clim <- a_clim * unique(env_mod[[3]][s]) # pick the Rsquared from this model and simulation
      
      elev_rank <- rank(sign(rho_elevation) * elev_z)
      clim_rank <- rank(env[,s]) # ranking the climate
      
      # NOTE: sign ensures ordering to be aligned with direction. low
      # ranks = low elevation, even with negative correaltion. this does not
      # matter of rho=positive, but for rh=negative
      score <- t(latent_ranks) +
        alpha * matrix(
          elev_rank,
          nrow = nnull,
          ncol = length(elev_rank),
          byrow = TRUE
        ) +
        alpha_clim * matrix(
          clim_rank,
          nrow = nnull,
          ncol = length(clim_rank),
          byrow = TRUE
        )
      
      score_ranks <- apply(score, 1, rank)
      occ_cf = apply(score_ranks, 2, function(x) occ_sorted[x])
      
      # test correlation
      occ_cf_elev_cor = apply(occ_cf, 2, function(x){cor(x, elev_z, method = "s")})
      
      ## variograms
      tmpvg <- data.table(x = xy[,1], y = xy[,2], occ_cf)
      sp::coordinates(tmpvg) <- ~x+y # turn into coordinates
      
      vg_sims <- data.table(
        iteration  = seq_len(ncol(occ_cf)),
        warning    = FALSE,
        range_param = NA_real_,
        cor_elevation = occ_cf_elev_cor
      )
      
      # check elevation correlation criterion FIRST (cheap) so we can skip
      # the expensive gstat fit for candidates that would be rejected anyway
      vg_sims[, keep_elev:=abs(cor_elevation - rho_elevation) < keep_tol_elev]
      vg_sims[, keep_sac:=FALSE]  # default; only set TRUE where we actually fit
      
      still_needed <- target_n - n_keep  # how many more acceptable columns needed
      
      for(vgi in seq_len(ncol(occ_cf))){
        
        # skip the variogram fit entirely if elevation criterion already fails
        if(!vg_sims$keep_elev[vgi]) next
        
        subtmp <- tmpvg[,vgi]
        names(subtmp) <- "response"
        vg <- gstat::variogram(response ~ 1, subtmp)
        fit <- withCallingHandlers(
          gstat::fit.variogram(vg, model = gstat::vgm(model = "Exp")),
          warning = function(w){
            vg_sims$warning[vgi] <<- TRUE
            invokeRestart("muffleWarning")
          }
        )
        vg_sims$range_param[vgi] <- fit$range[[1]]
        vg_sims$keep_sac[vgi] <- abs(vg_sims$range_param[vgi] - range_obs) / range_obs < keep_tol
        
        if(isTRUE(vg_sims$keep_sac[vgi])){
          still_needed <- still_needed - 1
          if(still_needed <= 0) break  # got enough from this batch, stop fitting
        }
      }
      
      # check which ones to keep, for range param and correlation with elevation
      keep_idx = which(vg_sims$keep_sac & vg_sims$keep_elev)
      
      # check if enough simulations meet requirements
      # check indexing!
      if(length(keep_idx) > 0){
        occ_keep <- cbind(occ_keep, occ_cf[, keep_idx, drop = FALSE])
        n_keep <- ncol(occ_keep)
      }
      batch <- batch + 1
    }
     
     ## store the 1 simulation (or NA if nothing can be found!)
     if(is.null(occ_keep)){
       cf[, s] = rep(NA, nrow(cf))
     }else{
       cf[, s] = occ_keep[,1:target_n, drop=FALSE] 
     }
     
    
     #cat("sim done:", "(", s, "/ 100", ")\n") 
    }
    
    saveRDS(cf, file = paste0("data/counterfactuals2/simulated_env_", species[i], "_", dataset_id, ".rds"))
    cat("species done:", "(", i, "/ ",length(species), ")\n")
    
  }, mc.cores = ncores)
}




