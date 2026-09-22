# counterfactual simulations analysis


# Load libraries ----------------------------------------------------------

# check libraries carfully. this is especially important on the server/cluster
.libPaths("/home/tietje/mamba/envs/r453/lib/R/library") 

rm(list=(ls())) # clear workspace
library(data.table)
library(parallel)
library(MASS)
library(fields)
library(sp)
library(RhpcBLASctl)



print(paste("System time: ", Sys.time()))


# Bash input ###############################

(args <- commandArgs(trailingOnly=TRUE))

dataset_ids   =as.numeric(strsplit(args[1], ",")[[1]])
which_CF      =args[2]
ncores        =as.numeric(args[3])
blas_thred    =as.numeric(args[4])
omp_thred     =as.numeric(args[5])


# example for local
dataset_ids   =17
which_CF      ="CF_IM"
ncores        =10
blas_thred    =1
omp_thred     =1


#################################################



# Pairwise analysis -------------------------------------------------------
# do this separately for each scenario, just to keep track of everything


for(dataset_id in dataset_ids){
 
   if(dataset_id==6){scalID=2}else{scalID=1}
  
  if(which_CF=="overlap"){

    files = dir("data/uncertainty", full.names = T, pattern = paste0(".*simulated_[A-Z].*_", dataset_id))
    nam = gsub(".*simulated_|_[0-9]{1,2}.*\\.rds", "", files)
    
    dat = readRDS(paste0("data/processed_occupancy_for_CFs_2_",dataset_id,".rds"))
    
    # creates pairs
    nam = nam[nam!=""] # no empties
    pairs = combn(nam, 2, simplify = T) # matrix object with ncol= n pairs and row=2 (each column one pair)
    gc()
    
    # restrict this to avoid BLAS thread oversubscription
    blas_set_num_threads(blas_thred)
    omp_set_num_threads(omp_thred)
    
    npairs = ncol(pairs)
    
    # create chunks for processing
    pair_chunks <- split(seq_len(npairs), cut(
      seq_len(npairs),
      breaks = ncores,
      labels = FALSE))
    
    # parallel process
    mclapply(seq_along(pair_chunks), function(chunk_id){
      res_list <- vector("list", length(pair_chunks[[chunk_id]]))
      for(i in seq_along(pair_chunks[[chunk_id]])){
        
        pair_id <- pair_chunks[[chunk_id]][i]
        sp <- pairs[, pair_id]
        tmp <- dat[scientificName %in% sp]
        
        ## check overlap
        tmp[, overlap_psi_03:=mean.psi[scientificName==sp[1]]>0.3 & mean.psi[scientificName==sp[2]]>0.3, by=.(siteID, endYear)]
        tmp[, overlap_psi_04:=mean.psi[scientificName==sp[1]]>0.4 & mean.psi[scientificName==sp[2]]>0.4, by=.(siteID, endYear)]
        tmp[, overlap_psi_05:=mean.psi[scientificName==sp[1]]>0.5 & mean.psi[scientificName==sp[2]]>0.5, by=.(siteID, endYear)]
        tmp[, overlap_psi_06:=mean.psi[scientificName==sp[1]]>0.6 & mean.psi[scientificName==sp[2]]>0.6, by=.(siteID, endYear)]
        tmp[, overlap_pa:=sum(pres.abs)>1, by=.(siteID, endYear)]
        tmp = tmp[, .(overlap_psi_03=any(overlap_psi_03),
                      overlap_psi_04=any(overlap_psi_04),
                      overlap_psi_05=any(overlap_psi_05),
                      overlap_psi_06=any(overlap_psi_06),
                      overlap_pa=any(overlap_pa)), .(endYear)]
        tmp$species_pair = paste0(sp, collapse = "|")
        res_list[[i]] = tmp
      }
      saveRDS(rbindlist(res_list), file = paste0("data/",dataset_id,"_overlap", "_chunk_", chunk_id, ".rds"))
      
    }, mc.cores = 5)
  }
  
  
  
  
 
  if(which_CF=="uncertainty"){
    
    files = dir("data/uncertainty", full.names = T, pattern = paste0(".*simulated_[A-Z].*_", dataset_id))
    nam = gsub(".*simulated_|_[0-9]{1,2}.*\\.rds", "", files)
  
    dat = readRDS(paste0("data/processed_occupancy_for_CFs_2_",dataset_id,".rds"))
    
    # creates pairs
    nam = nam[nam!=""] # no empties
    pairs = combn(nam, 2, simplify = T) # matrix object with ncol= n pairs and row=2 (each column one pair)
    gc()
    
    # restrict this to avoid BLAS thread oversubscription
    blas_set_num_threads(blas_thred)
    omp_set_num_threads(omp_thred)
    
    npairs = ncol(pairs)
    
    # create chunks for processing
    pair_chunks <- split(seq_len(npairs), cut(
      seq_len(npairs),
      breaks = ncores,
      labels = FALSE))
    
    # parallel process
    mclapply(seq_along(pair_chunks), function(chunk_id){
      res_list <- vector("list", length(pair_chunks[[chunk_id]]))
      for(i in seq_along(pair_chunks[[chunk_id]])){
        
        pair_id <- pair_chunks[[chunk_id]][i]
        sp <- pairs[, pair_id]
        tmp <- dat[scientificName %in% sp]
        
        ## observed correlation
        cor_org = tmp[, .(rho = cor.test(mean.psi[scientificName==sp[1]], mean.psi[scientificName==sp[2]], 
                                         method="s")$estimate,
                          rho_pvalue = cor.test(mean.psi[scientificName==sp[1]], mean.psi[scientificName==sp[2]], 
                                                method="s")$p.value,
                          #c_score = bipartite::C.score(cbind(pres.abs[scientificName==pairs[2, i]], pres.abs[scientificName==pairs[2, i]]), normalise=TRUE),
                          iteration = 1,
                          treatment = "original",
                          n = 1
                          
        ),
        by=.(endYear)]
        cor_org$species_pair = paste(sp[1], sp[2], sep="|")
        
        ## load simulations (basename is an exact match)
        u_list1 <- readRDS(files[grep(sp[1], files)])
        u_list2 <- readRDS(files[grep(sp[2], files)])
        
        ## correlations from uncertainty distributions
        ## T1
          cor_u1 <- rbindlist(lapply(seq_len(ncol(u_list1[[1]])), function(m){
            
            ct <- cor.test(u_list1[[1]][,m], u_list2[[1]][,m], method = "s")
            data.table(rho = unname(ct$estimate), 
                       rho_pvalue = ct$p.value,
                       iteration = m,
                       endYear = 1,
                       treatment = which_CF,
                       n = 100)
          }))
          cor_u2 <- rbindlist(lapply(seq_len(ncol(u_list1[[2]])), function(m){
            
            ct <- cor.test(u_list1[[2]][,m], u_list2[[2]][,m], method = "s")
            data.table(rho = unname(ct$estimate), 
                       rho_pvalue = ct$p.value,
                       iteration = m,
                       endYear = 3,
                       treatment = which_CF,
                       n = 100)
          }))
        cor_u <- rbind(cor_u1, cor_u2)
        cor_u[, species_pair:=paste(sp[1], sp[2], sep="|")]
        
        # rbind results original T1 + T2 and simulations
        res <- rbind(cor_org, cor_u)
        res[, dataset_id:=dataset_id]
        
        # safe in list
        res_list[[i]] <- res
    }
      
      ## save one file per chunk
      saveRDS(rbindlist(res_list), 
              file = paste0("data/uncertainty/results/results_",dataset_id,"_",which_CF,"_chunk_", chunk_id, ".rds"))
      NULL
      
    }, mc.cores = ncores)
  }  
  
  
  
  
  dat = readRDS(paste0("data/processed_occupancy_for_CFs_2_",dataset_id,".rds"))
  
  if(which_CF=="CF_IM"){
    files = dir("data/counterfactuals2", full.names = T, pattern = paste0(".*simulated_[A-Z].*_", dataset_id))
    files = files[grep(paste0(unique(dat$scientificName), collapse = "|"), files)]
    nam = gsub(".*simulated_|_[0-9]{1,2}.*\\.rds", "", files)
  } 
  if(which_CF=="CF_IM_elev"){
    files = dir("data/counterfactuals2", full.names = T, pattern = paste0(".*simulated_elevation_.*_", dataset_id))
    files = files[grep(paste0(unique(dat$scientificName), collapse = "|"), files)]
    nam = gsub(".*simulated_elevation_|_[0-9]{1,2}.*\\.rds", "", files)
  } 
  if(which_CF=="CF_IM_env"){
    files = dir("data/counterfactuals2", full.names = T, pattern = paste0(".*simulated_env.*_", dataset_id))
    files = files[grep(paste0(unique(dat$scientificName), collapse = "|"), files)]
    nam = gsub(".*simulated_env_|_[0-9]{1,2}.*\\.rds", "", files)
  } 
  if(which_CF=="clim"){
    files = dir("data/counterfactuals2", full.names = T, pattern = paste0(".*c\\+l.*_", dataset_id))
    files = files[grep(paste0(unique(dat$scientificName), collapse = "|"), files)]
    nam = gsub(".*species_|_[0-9]{1,2}.*\\.rds", "", files)
  } 
  
  
  # creates pairs
  nam = nam[nam!=""] # no empties
  pairs = combn(nam, 2, simplify = T) # matrix object with ncol= n pairs and row=2 (each column one pair)
  gc()
  
  # restrict this to avoid BLAS thread oversubscription
  blas_set_num_threads(blas_thred)
  omp_set_num_threads(omp_thred)
  
  npairs = ncol(pairs)

  # create chunks for processing
  pair_chunks <- split(seq_len(npairs), cut(
    seq_len(npairs),
    breaks = ncores,
    labels = FALSE))

  # parallel process
  mclapply(seq_along(pair_chunks), function(chunk_id){
    res_list <- vector("list", length(pair_chunks[[chunk_id]]))
    for(i in seq_along(pair_chunks[[chunk_id]])){
      
      pair_id <- pair_chunks[[chunk_id]][i]
      sp <- pairs[, pair_id]
      tmp <- dat[scientificName %in% sp]
      
      ## observed correlation
      cor_org = tmp[, .(rho = cor.test(mean.psi[scientificName==sp[1]], mean.psi[scientificName==sp[2]], 
                                       method="s")$estimate,
                        rho_pvalue = cor.test(mean.psi[scientificName==sp[1]], mean.psi[scientificName==sp[2]], 
                                              method="s")$p.value,
                        #c_score = bipartite::C.score(cbind(pres.abs[scientificName==pairs[2, i]], pres.abs[scientificName==pairs[2, i]]), normalise=TRUE),
                        iteration = 1,
                        treatment = "original",
                        n = 1
      ),
      by=.(endYear)]
      cor_org$species_pair = paste(sp[1], sp[2], sep="|")
      
      if(which_CF=="clim"){
        cf_list1 <- readRDS(files[grep(sp[1], files)])
        cf_list2 <- readRDS(files[grep(sp[2], files)])
        cf_list1 = cf_list1[[1]]
        cf_list2 = cf_list2[[1]]
        
      }else{
        cf_list1 <- readRDS(files[grep(sp[1], files)])
        cf_list2 <- readRDS(files[grep(sp[2], files)])
      }
      ## load simulations (basename is an exact match)

      
      ## simulated correlations
      # this is a bit tricky now since we need to sample for partner of a
      # species that has < 100 simulations
      shorty = min(ncol(cf_list1), ncol(cf_list2))
       # species that has < 100 simulations.
      # some even have zero
      shorty = min(ncol(cf_list1), ncol(cf_list2))
      
      if(is.null(ncol(cf_list1))|is.null(ncol(cf_list2))){
        cor_cf = data.table(rho = NA, 
                            rho_pvalue = NA,
                            iteration = 1,
                            endYear = 3,
                            treatment = which_CF,
                            n = 0)
      }else{
        cor_cf <- rbindlist(lapply(seq_len(shorty), function(m){
          
          ct <- cor.test(cf_list1[,m], cf_list2[,m], method = "s")
          data.table(rho = unname(ct$estimate), 
                     rho_pvalue = ct$p.value,
                     iteration = m,
                     endYear = 3,
                     treatment = which_CF,
                     n = shorty)
        }))
      }
      cor_cf[, species_pair:=paste(sp[1], sp[2], sep="|")]
      
      # rbind results original T1 + T2 and simulations
      res <- rbind(cor_org, cor_cf)
      res[, dataset_id:=dataset_id]
      
      # safe in list
      res_list[[i]] <- res
    }
    print(i)
    ## save one file per chunk
    saveRDS(rbindlist(res_list), 
            file = paste0("data/counterfactuals2/results/results_",dataset_id,"_",which_CF,"_chunk_", chunk_id, ".rds"))
    NULL
    
  }, mc.cores = ncores)
  
}






