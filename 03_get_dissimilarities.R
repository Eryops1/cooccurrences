# Script: 03_get_dissimilarities.R
# Author: Melanie Tietje
# Email: tietje@fzp.czu.cz
# GitHub: @Eryops1
# Last Modified: 2025-10-30
# Purpose: Calculate dissimilarities in co-occurrence measures between sampling
#          periods, calculate mantel tests for influence of phylogeny and traits.
# Output: .rds objects named "data/processed_spass.rds",
#          "sum_mean_psi_for_maps.rds", "range_change.rds", "mantel_results.rds",
#          "braycurtis.rds"
# Notes: groundhog will ensure the exact same R packages will be used and create
#        a library on first run, which might take a while

# Load libraries ----------------------------------------------------------

rm(list=(ls())) # clear workspace
library(groundhog)
pkgs <- c("data.table",
          "sf",
          "clootl",
          "ape",
          "ecodist",
          "vegan")
groundhog.library(pkgs, "2025-10-25")







# read data --------------------------------------------------------------------


## occupancy data ----- 
dataset_id =c(5,6,17,26)
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
  tmp[[i]]$datasetID = gsub(".*occupancy_|_.*\\.rds", "", files[i])
  tmp[[i]]$scalingID = 1
  tmp[[i]]$endYear = gsub(".*_|\\.rds", "", files[i])
}

dat = rbindlist(tmp)


# set scale to 2 for NY
if(6 %in% dataset_id){
  dat$scalingID[dat$datasetID==6]=2
}

# subset to species that occur in all time periods:
## step 1: remove all NAs
dat = na.omit(dat)
## step2: count years with data
dat[, alltime:=length(unique(endYear)), by=.(datasetID, scientificName)]
table(dat$alltime, dat$datasetID)

# delete species that do not meet required number of sample periods
cond = as.data.table(table(gsub(".*occupancy_|_.*\\.rds", "", files)))
dat = dat[datasetID==cond$V1[1] & alltime==cond$N[1] | 
            datasetID==cond$V1[2] & alltime==cond$N[2] | 
            datasetID==cond$V1[3] & alltime==cond$N[3] | 
            datasetID==cond$V1[4] & alltime==cond$N[4]]
table(dat$alltime, dat$datasetID)
any(is.na(dat$mean.psi)) # no NAs? good!

# save sum occupancy probabilities for making maps later
dat[, sum.psi_siteID:=sum(mean.psi),by=.(siteID, datasetID, endYear)]
maps = unique(dat[,.(sum.psi_siteID, siteID, datasetID, endYear)])
saveRDS(maps, "data/sum_mean_psi_for_maps.rds")

# get changes in ranges 
# change = sum occupancy last time - sum occupancy first time
dat[, mean.psi_sum:=sum(mean.psi), by=.(datasetID, scientificName, endYear)]
dat[, change_occupancy:=(sum(mean.psi[endYear==max(endYear)]) - sum(mean.psi[endYear==1])) /sum(mean.psi[endYear==1]), 
    by=.(datasetID, scientificName)]

# get change in present grid cells
dat[, pres.abs_sum:=sum(pres.abs), by=.(datasetID, scientificName, endYear)]
dat[, change_pres.abs:=(sum(pres.abs[endYear==max(endYear)]) - sum(pres.abs[endYear==1])) /sum(pres.abs[endYear==1]), 
    by=.(datasetID, scientificName)]

dat$atlas = dat$datasetID
dat$atlas = gsub("26", "Europe", dat$atlas)
dat$atlas = gsub("5", "Czechia", dat$atlas)
dat$atlas = gsub("6", "New York", dat$atlas)
dat$atlas = gsub("17", "New Zealand", dat$atlas)

saveRDS(unique(dat[, .(datasetID, endYear, scientificName, scalingID, alltime, 
                       mean.psi_sum, pres.abs_sum, change_occupancy, change_pres.abs)]), 
               "data/range_change.rds")





## phylo data -----
phy = clootl::extractTree()

any(phy$edge.length==0)
ape::is.binary(phy)
ape::is.rooted(phy)
ape::is.ultrametric(phy, tol = .Machine$double.eps^0.4) # ok




## avonet -----
avonet = fread("data/AVONET1_BirdLife.csv")
avonet$species = gsub(" ", "_", avonet$Species1)

# select traits to use
input <- avonet[,.(species, Beak.Length_Culmen, Beak.Width, Beak.Depth, Tarsus.Length, Wing.Length, 
                   Kipps.Distance, Secondary1, `Hand-Wing.Index`, Tail.Length, Mass)]




## Match data -----------------------------------------------------------------
atlas = sort(unique(dat$datasetID))  
comm_sp = list()
for(a in 1:length(atlas)){
  p = intersect(unique(dat$scientificName[dat$datasetID==atlas[a]]), phy$tip.label)
  f = intersect(unique(dat$scientificName[dat$datasetID==atlas[a]]), input$species)
  comm_sp[[a]] =  intersect(p,f)
}
names(comm_sp) = atlas

dat = dat[scientificName %in% comm_sp$`5` & datasetID=="5" |
                scientificName %in% comm_sp$`6` & datasetID=="6"| 
            scientificName %in% comm_sp$`26` & datasetID=="26" |
          scientificName %in% comm_sp$`17` & datasetID=="17", ] 
gc()













# Get correlation matrix -------------------------------------------------------
# build Spearman correlation matrix from occupancy data for Mantel tests

mat = list()
atlas_vec = sort(unique(dat$datasetID))
for(a in 1:length(atlas_vec)){
  tmp = dat[datasetID %in% atlas_vec[a],]
  times = sort(unique(tmp$endYear))
  mat_int = list()
  for(t in 1:length(times)){
    tmp2 = tmp[endYear==times[t],]
    tmp2 = dcast(tmp2, rowid(scientificName) + siteID ~ scientificName, value.var = 'mean.psi', fill = NA)
    tmp2 = tmp2[, -c("scientificName")]
    tmp2 = as.dist(cor(tmp2[, -c("siteID")], method="s"))
    mat_int[[t]] = list(data.table(dataset_id = atlas_vec[a],
                                   time_bin = t),
                        tmp2)
  }
  mat[[a]] = mat_int
}
names(mat) = atlas

# get difference in correlation between sampling periods
dd = list()
for(i in 1:length(mat)){
  tmp = mat[[i]]
  dd[[i]] = abs(tmp[[length(tmp)]][[2]] - tmp[[1]][[2]])
}
names(dd) = atlas

# saveRDS(dd, "data/distance_matrix_abs.rds")
dd = readRDS("data/distance_matrix_abs.rds")



## PD matrix ------------------------------------------------------------------- 
pd = list()
for(a in 1:length(atlas)){
  sp = sort(unique(dat$scientificName[dat$datasetID==atlas[a]]))
  keep = sp[sp %in% phy$tip.label]
  phy.sub <- ape::keep.tip(phy, keep)
  
  tmp = ape::cophenetic.phylo(phy.sub)
  pd[[a]] = as.dist(tmp)
}
names(pd) = atlas

#saveRDS(pd, "data/PD_distance_matrix.rds")



## FD matrix ------------------------------------------------------------------
fd = list()
for(a in 1:length(atlas)){
  sp = sort(unique(dat$scientificName[dat$datasetID==atlas[a]]))
  keep = sp[sp %in% input$species]
  input.sub <- input[species %in% keep,]
  
  any(is.na(input))
  input.sub = na.omit(input.sub)
  gow = cluster::daisy(input.sub[,!"species"], metric = "gower")
  
  gow = as.matrix(gow)
  colnames(gow) = input.sub$species
  row.names(gow) = colnames(gow)
  
  fd[[a]] = as.dist(gow)
}
names(fd) = atlas

#saveRDS(fd, "data/FD_distance_matrix.rds")




# Mantel test ------------------------------------------------------------------

## PD -------------
if(any(names(dd)!=names(pd))){print("STOP!! order not correct")}

atli <- names(dd)
mantel_PD_res <- list()
mrm_result = list()
for(a in 1:length(atli)){
  sub_dd <- dd[[a]]
  sub_pd <- pd[[a]]
  # rearrange PD distance matrix to match order in dd
  dist_matrix1 <- as.matrix(sub_dd)
  dist_matrix <- as.matrix(sub_pd)
  # reorder
  new_order <- match(colnames(dist_matrix1), colnames(dist_matrix))
  reordered_matrix <- dist_matrix[new_order, new_order]
  # back to dist object
  pd_reord =  as.dist(reordered_matrix)
  
  names(sub_dd) == names(pd_reord)
  rownames(sub_dd) == rownames(pd_reord)
  # nice!
  
  mantel_PD_res[[a]] = vegan::mantel(sub_dd, pd_reord, method = "s")
  
  # alt: MRM (multiple regression on distance matrices)
  # mrm_result[[a]] = MRM(sub_dd ~ pd_reord, nperm = 999)
}

names(mantel_PD_res) <- atli
#names(mrm_result) <- atli


mantel_PD_res$`5`
#mrm_result$`5`






## FD -------------
if(any(names(dd)!=names(fd))){print("STOP!! order not correct")}

atli <- names(dd)
mantel_FD_res <- list()
for(a in 1:length(atli)){
  sub_dd <- dd[[a]]
  sub_fd <- fd[[a]]
  # rearrange PD distance matrix to match order in dd
  dist_matrix1 <- as.matrix(sub_dd)
  dist_matrix <- as.matrix(sub_fd)
  # reorder
  new_order <- match(colnames(dist_matrix1), colnames(dist_matrix))
  reordered_matrix <- dist_matrix[new_order, new_order]
  # back to dist object
  fd_reord =  as.dist(reordered_matrix)
  
  names(sub_dd) == names(fd_reord)
  rownames(sub_dd) == rownames(fd_reord)
  # nice!
  
  mantel_FD_res[[a]] = vegan::mantel(sub_dd, fd_reord, method = "s")
}

names(mantel_FD_res) <- atli



## save ----
tmp = list(mantel_PD_res, mantel_FD_res)
names(tmp) = c('mantel_PD_res', 'mantel_FD_res')
saveRDS(tmp, "data/mantel_results.rds")













# Jaccard and BrayCurtis for species-species -----------------------------------

# Get Bray-Curtis for Fig 3c (see also Fig S3)

# Step 1: get BrayCurtis matrices (this is similar to Spearman correlation, but
# avoids negative values)

atlas_list = list()
atli = sort(unique(dat$datasetID))
mat = list()
for(a in 1:length(atli)){
  tmp = dat[datasetID==atli[a],]
  times = sort(unique(tmp$endYear))
  mat_int = list()
  for(t in 1:length(times)){
    tmp2 = tmp[endYear==times[t],]
    tmp2 = dcast(tmp2, rowid(scientificName) + siteID ~ scientificName, value.var = 'mean.psi', fill = NA)
    tmp2 = tmp2[, -c("scientificName")]
    tmp2 = vegdist(t(tmp2[, -c("siteID")]), method="bray", upper=TRUE, diag=TRUE)
    mat_int[[t]] = list(data.table(dataset_id = atli[a],
                                   time_bin = t), 
                        tmp2)
  }
  mat[[a]] = mat_int
}
names(mat) = atli


# Step 2: get BrayCurtis comparing BrayCurtis dissimilarity with other species
# across sampling periods, for each focal species. Think about it as comparing
# the composition of spatial overlaps that a species has.

bc_res = list()
for(a in 1:length(atli)){
  time1 = mat[[a]][[1]][[2]]
  time2 = mat[[a]][[length(mat[[a]])]][[2]]
  print(all(names(time1)==names(time2))) # test if names match
  tmp = data.table(species=names(time1),
                   bray_curtis = NA, 
                   atlas=atli[a])
  time1_m = as.matrix(time1)
  time2_m = as.matrix(time2)
  for(l in 1:ncol(time1_m)){
    bc = vegdist(t(data.frame(time1_m[,l], time2_m[,l])))
    tmp$bray_curtis[l] = as.numeric(bc)
  }
  bc_res[[a]] = tmp
}

bc = rbindlist(bc_res)
saveRDS(bc, "data/braycurtis.rds")

















#  # OLD JACCARD VERSION
# jaccard <- function(a, b) {
#   intersection = length(intersect(a, b))
#   union = length(a) + length(b) - intersection
#   return (intersection/union)
# }
# 
# 
# atlas_list = list()
# atli = sort(unique(dat$datasetID))
# for(a in 1:length(atli)){         # for each atlas
#   tmp_a = dat[datasetID==atli[a],]
#   times = sort(unique(tmp_a$endYear))
#   times_list = list()
#   for(t in 1:length(times)){                      # for each time bin
#     tmp_t = tmp_a[endYear==times[t],]
#     mat <- dcast(tmp_t, siteID ~ scientificName, value.var = "pres.abs")
#     mat = mat[,-"siteID"]
#     sp = names(mat)
#     coocc <- list()
#     for(i in 1:length(sp)){
#       target_sp = sp[i]
#       # reduce to rows where target sp exists
#       tmp <- mat[which(mat[,get(target_sp)]==1),]
#       # get species names whichs col sums are not 0
#       friends <- names(which(colSums(tmp)>0))
#       # remove target species name
#       friends = friends[!grepl(target_sp, friends)]
#       coocc[[i]] = friends
#     }
#     names(coocc) = sp
#     times_list[[t]] = coocc
#   }
#   names(times_list) = times
#   
#   # get jaccard for each species between first and last time bin
#   sp = names(times_list[[1]])
#   jac = data.table(species = sp,
#                    jaccard = NA, 
#                    friends = NA)
#   for(j in 1:length(sp)){
#     friends1 = times_list[[times[1]]][[sp[j]]]
#     friends2 = times_list[[tail(times,1)]][[sp[j]]]
#     jac$jaccard[j] = jaccard(friends1, friends2)
#     jac$scientificName[j] = sp[j]
#     jac$friendsT1[j] = length(friends1)
#     jac$friendsTend[j] = length(friends2)
#   }
#   # attach atlas ID
#   jac$datasetID = atli[a]
#   # save everything
#   atlas_list[[a]] = jac 
#   cat(a, "\r")
# }
# 
# atlas_list
# jac = rbindlist(atlas_list)
# saveRDS(jac, "data/jaccard.rds")
# 
# 
# 




