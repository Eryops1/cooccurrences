# Main analysis after the computational counterfactual stages.
#
# Required upstream outputs:
#   - counterfactual pairwise results
#   - uncertainty results and overlap files
#   - data/z_scores.rds from 03_z_scores_server.R
#   - data/trait_analysis_input2.rds

rm(list = ls())

library(data.table)
library(ggplot2)
library(terra)
library(sf)
library(cowplot)
library(gstat)
library(sp)
library(ggExtra)
library(ggpubr)
library(ggpmisc)
library(ggsci)
library(scico)
library(circlize)
library(mgcv)
library(clootl)
library(ape)
library(cluster)
library(phytools)
library(performance)
library(lmerMultiMember)
library(knitr)
library(broom)
library(broom.mixed)

theme_set(theme_classic(base_size = 9)+
            theme(strip.background = element_blank()))
            # theme(axis.title = ggtext::element_markdown())+
            # theme(axis.title.x = ggtext::element_markdown())+
            # theme(axis.title.y = ggtext::element_markdown())+
            # theme(legend.title = ggtext::element_markdown()))



source("99_functions.R")





#----------------------------------------------------------------------------#
# load data -----------------------------------------------------------------
#----------------------------------------------------------------------------#


files = dir("data/counterfactuals2/results", full.names = T, pattern = "results_.*chunk")
res = lapply(files, readRDS)
res = rbindlist(res, fill=TRUE, use.names = TRUE)
table(res$treatment)

res$treatment = factor(res$treatment, levels=c("original", "clim", "CF_IM", "CF_IM_elev", "CF_IM_env"))

# replace larger end Year with 2 for consistency
res[endYear==max(endYear), "endYear"] = 2


# check some basics, remove dups since the original correlation is in there per treatment
names(res)
dup = duplicated(res, by=c("endYear", "rho", "rho_pvalue", "iteration",
                           "treatment", "species_pair", "dataset_id")) 
res = res[-which(dup),]
rm(dup)

# also get rid of pairs that just have no simulations
res <- res[n!=0, ]

# this is gonna take a while, but we gotta make sure all is in order 
#apply(res, 2, function(x){any(is.na(x))})

(rem = res[which(is.na(res$rho)), ]) # remove species pairs if missing rho anywhere
res = res[!(species_pair %in% rem$species_pair & dataset_id%in%rem$dataset_id),]


# only keep pairs that are present in all scenarios
tmp <- res[, .N, .(species_pair, dataset_id)] # should be 402, anything less --> out
table(tmp$N)
# remove all pairs with less than 402
rem <- tmp$species_pair[which(tmp$N!=402)]
res <- res[!species_pair %in% rem, ]





# uncertainty
files = dir("data/uncertainty/results", full.names = T, pattern = "^results_.*chunk")
un = lapply(files, readRDS)
un = rbindlist(un, fill=TRUE, use.names = TRUE)
un$treatment = factor(un$treatment, levels=c("original", "uncertainty"))

# subset to the pairs present in res
pair_at_comb = unique(paste0(res$species_pair, res$dataset_id))
un = un[paste0(species_pair, dataset_id) %in% pair_at_comb, ]
un = un[treatment=="uncertainty",]
un[endYear==max(endYear), "endYear"] = 2



# merge into res
res <- rbind(res, un)
res = unique(res)

gc()


# Get delta rho 

## collect rho_t1 per treatment
res[treatment=="original", rho_t1:=rho[endYear==1], by=.(species_pair, dataset_id)]
res[treatment!="original", rho_t1:=rho[endYear==1 & treatment=="uncertainty"], by=.(species_pair, dataset_id, iteration)]
any(is.na(res$rho_t1))

res[, delta_cor:=rho - rho_t1, .(species_pair, dataset_id, treatment, iteration)]


# remove the non overlapping ones - this is based on the average species maps
files = dir("data", pattern = "overlap", full.names = T)
over = lapply(files, readRDS)
names(over) = gsub("data/|_.*", "", files)
over = rbindlist(over, idcol = "dataset_id")
over[endYear==3, 'endYear'] = 2

over[, overlap_pa_type:=ifelse(
  overlap_pa[endYear==1] & !overlap_pa[endYear==2], "only_t1", ifelse(
    !overlap_pa[endYear==1] & overlap_pa[endYear==2], "only_t2", ifelse(
      overlap_pa[endYear==1] & overlap_pa[endYear==2], "both", "none")
  )), .(species_pair, dataset_id)]

cols <- grep("^overlap_psi_0[3-6]$", names(over), value = TRUE)
for (col_in in cols) {
  col_out <- paste0(col_in, "_type")
  over[, (col_out) := {
    v <- get(col_in)
    fifelse(
      v[endYear == 1] & !v[endYear == 2], "only_t1",
      fifelse(
        !v[endYear == 1] & v[endYear == 2], "only_t2",
        fifelse(v[endYear == 1] & v[endYear == 2], "both", "none")
      ))}, by = .(species_pair, dataset_id)]
}

keep = unique(over[, .(species_pair, dataset_id, overlap_psi_03_type)])
keep$dataset_id = as.numeric(keep$dataset_id)
res = merge(res, keep, all.x=T)
res = res[overlap_psi_03_type %in% c("both", "only_t1", "only_t2"),]



# reduce to the last year, T1 rho is stored in column, we dont need T1 rows anymore
saveRDS(res[treatment=="uncertainty"], "output/corT1T2.rds")
table(res$endYear)
res <- res[endYear==max(endYear),]



knitr::kable(res[treatment=="original", .N, .(dataset_id, endYear)])

res[, c("sp1", "sp2") := tstrsplit(species_pair, "|", fixed = TRUE)]
spcount <- unique(res[, .(sp1,sp2,treatment,dataset_id)])
spcount <- melt(spcount, measure.vars = c("sp1","sp2"))
spcount <- unique(spcount[, .(treatment, dataset_id, value)])
knitr::kable(spcount[, .N, .(dataset_id, treatment)])


res$treatment = as.character(res$treatment)
res$treatment[res$treatment=="clim"] = "clim_lu"
res$treatment[res$treatment=="original"] = "average_model"
res$treatment = factor(res$treatment, 
                       levels = c("average_model","clim_lu","CF_IM","CF_IM_elev","CF_IM_env","uncertainty"),
                       labels = c("AM", "S1", "S2", "S3", "S4", "observed"))


## add stats for distribution differences from zero
# p value great = the proportion of values larger or equal to zero in the total
# proportions. similar for other site
res[, `:=`(
  p.value_greater = (sum(delta_cor <= 0) + 1) / (.N + 1),
  p.value_less    = (sum(delta_cor >= 0) + 1) / (.N + 1)
  
), .(dataset_id, species_pair, treatment)]
res[, p_two_sided := pmin(2 * pmin(p.value_greater, p.value_less), 1)]





# pairwise values
pairwise_res = unique(res[, .(delta_cor_median=median(delta_cor),
                       delta_cor_mad=mad(delta_cor),
                       ymin = min(delta_cor),
                       lower025 = quantile(delta_cor, 0.025),
                       lower = quantile(delta_cor, 0.25),
                       middle = median(delta_cor),
                       upper = quantile(delta_cor, 0.75),
                       upper975 = quantile(delta_cor, 0.975),
                       ymax = max(delta_cor),
                       p.greater_zero = p.value_greater,
                       p.less_zero = p.value_less
), 
.(treatment, species_pair, dataset_id)])

pairwise_res[, delta_cor_obs:=delta_cor_median[treatment=="observed"], .(dataset_id, species_pair)]
pairwise_res[, greater_less:=ifelse(p.greater_zero<0.05, "greater", ifelse(p.less_zero<0.05, "less", "neutral"))]

pairwise_res[, y_key:=rank(middle), .(treatment, dataset_id)] # key for caterpillar plot y axis




saveRDS(list(res, pairwise_res), "data/final_analysis_input2.rds")














#----------------------------------------------------------------------------#
# Qualitative association transitions
#----------------------------------------------------------------------------#

zscores <- readRDS("data/z_scores.rds")
zscores <- rbindlist(zscores)
zscores[, iteration:=rep(c(1:100)), .(dataset_id, timestep, species_pair)]
zscores$atlas = factor(zscores$dataset_id, levels = c(5,6,17,26),
                   labels=c("Czechia", "New York", "New Zealand", "Europe"))
zscores$timestep[zscores$timestep==3] = 2

zscores[, state:=ifelse(z_scores<(-1.96), "segregated", ifelse(z_scores>1.96, "aggregated", "neutral")), 
      .(species_pair, timestep, atlas, iteration)]

z_wide = dcast(zscores, species_pair + atlas + iteration ~ timestep, value.var = "state")
z_wide[, transition:=paste(`1`, `2`, sep=">>")]

circ_in <- z_wide[, .N, .(atlas, transition, iteration)][, perc := 100 * N / sum(N),  .(atlas, iteration)]
# total % leaving / staying in each category = sum of perc where category is the source
circ_in[, from:=gsub(">>.*", "", transition)]
circ_in[, to:=gsub(".*>>", "", transition)]

cat_totals_from = unique(circ_in[, .(total = sum(perc)), .(atlas, from, iteration)])
cat_totals_from = cat_totals_from[, .(median = median(total, na.rm=T),
                    mad = mad(total, na.rm=T),
                    q025=quantile(total, 0.025),
                    q975=quantile(total, 0.975)), .(atlas, from)]
cat_totals_remain = unique(circ_in[, .(total = sum(perc[from==to])), .(atlas, from, iteration)])
cat_totals_remain = cat_totals_remain[, .(median = median(total, na.rm=T),
                                      mad = mad(total, na.rm=T),
                                      q025=quantile(total, 0.025, na.rm=T),
                                      q975=quantile(total, 0.975, na.rm=T)), .(atlas, from)]
cat_totals_leave = unique(circ_in[, .(total = sum(perc[from!=to])), .(atlas, from, iteration)])
cat_totals_leave = cat_totals_leave[, .(median = median(total, na.rm=T),
                                          mad = mad(total, na.rm=T),
                                          q025=quantile(total, 0.025, na.rm=T),
                                          q975=quantile(total, 0.975, na.rm=T)), .(atlas, from)]


# add summary across iterations
circ_in_avg <- unique(circ_in[, .(
  median_N = median(N, na.rm=T),
  median_perc = median(perc, na.rm=T),
  mad_perc = mad(perc, na.rm=T),
  q025=quantile(perc, 0.025),
  q975=quantile(perc, 0.975),
  trans=ifelse(from==to, "stable", "changing")), .(atlas, transition)])
circ_in_avg$transition <- factor(circ_in_avg$transition, levels = c("aggregated>>aggregated",
                                                          "segregated>>segregated",
                                                          "neutral>>neutral",
                                                          
                                                          "neutral>>aggregated",
                                                          "neutral>>segregated",
                                                          "aggregated>>neutral",
                                                          "segregated>>neutral",
                                                          
                                                          "aggregated>>segregated",
                                                          "segregated>>aggregated"
))

# stats --
circ_in_tmp = circ_in_avg[, .(perc_N=sum(median_N),
                              perc_trans=sum(median_perc)), .(atlas, trans)]
knitr::kable(circ_in_tmp[order(atlas),.(atlas, trans, round(perc_N), perc_trans),], digits = 2)

circ_in_tmp[, .(mean_perc=mean(perc_trans)), .(atlas, trans)]
circ_in_tmp[, .(mean_perc=mean(perc_trans)), .(trans)]

# # numbers, not percent, based on all iterations
# circ_in_tmp = circ_in[, .(N_trans=sum(N)), .(atlas, trans)]
# tstat = circ_in_tmp[, .(sum_N=sum(N_trans)), .(trans)]
# tstat$sum_N[1]/sum(tstat$sum_N)

knitr::kable(circ_in_avg[order(atlas, trans),.(atlas, transition, trans, median_N, median_perc, mad_perc)], digits = 3)
fwrite(circ_in_avg, "output/transition_stats_last.csv")

# Save qualitative-transition analysis products.
saveRDS(list(zscores = zscores, z_wide = z_wide, circ_in = circ_in,
            cat_totals_from = cat_totals_from,
            cat_totals_remain = cat_totals_remain,
            cat_totals_leave = cat_totals_leave,
            circ_in_avg = circ_in_avg),
        "data/qualitative_transition_analysis.rds")
















#----------------------------------------------------------------------------#
# Trait analysis
#----------------------------------------------------------------------------#


res = readRDS("data/trait_analysis_input2.rds")
res = res[treatment %in% c('uncertainty'), 
          .(dataset_id, species_pair, treatment, delta_cor, iteration)]

res[, delta_cor_median:=median(delta_cor), .(species_pair, dataset_id, treatment)]
res[, delta_cor_mad:=mad(delta_cor), .(species_pair, dataset_id, treatment)]
res = unique(res[, .(dataset_id, species_pair, treatment, delta_cor_median, delta_cor_mad)])



## get PD ------------------------------------------------------------------

phy = clootl::extractTree(species="all_species")
phy$tip.label = gsub(" ", "_", phy$tip.label)
any(phy$edge.length==0)
ape::is.binary(phy)
ape::is.rooted(phy)
ape::is.ultrametric(phy, tol = .Machine$double.eps^0.4) # ok


pairs = sort(unique(res$species_pair))
pair_dt <- tstrsplit(pairs, "|", fixed = TRUE)

sp1 <- pair_dt[[1]]
sp2 <- pair_dt[[2]]

# Keep only species actually needed
spp <- unique(c(sp1, sp2))

## tax matching avibase (28 Juli 2026) --------------------------
spp[which(!spp %in% phy$tip.label)]

# Accipiter gentilis,  	Astur gentilis in bird life
"Astur_gentilis" %in% phy$tip.label
phy$tip.label[which(phy$tip.label=="Astur_gentilis")] = "Accipiter_gentilis"

# Cahalcites_lucidus, spelling: Chalcites_lucidus
"Chalcites_lucidus" %in% phy$tip.label
phy$tip.label[which(phy$tip.label=="Chalcites_lucidus")] = "Cahalcites_lucidus"

# Charadrius_alexandrinus = Anarhynchus alexandrinus
"Anarhynchus_alexandrinus" %in% phy$tip.label
phy$tip.label[which(phy$tip.label=="Anarhynchus_alexandrinus")] = "Charadrius_alexandrinus"

# Charadrius_bicinctus = Anarhynchus_bicinctus
"Anarhynchus_bicinctus" %in% phy$tip.label
phy$tip.label[which(phy$tip.label=="Anarhynchus_bicinctus")] = "Charadrius_bicinctus"

#Charadrius_dubius = Thinornis dubius
"Thinornis_dubius" %in% phy$tip.label
phy$tip.label[which(phy$tip.label=="Thinornis_dubius")] = "Charadrius_dubius"

#Corvus_monedula = Coloeus_monedula 
"Coloeus_monedula" %in% phy$tip.label
phy$tip.label[which(phy$tip.label=="Coloeus_monedula")] = "Corvus_monedula"

# Ixobrychus_exilis =  	Botaurus_exilis
"Botaurus_exilis" %in% phy$tip.label
phy$tip.label[which(phy$tip.label=="Botaurus_exilis")] = "Ixobrychus_exilis"

# Ixobrychus_minutus = Botaurus_minutus
"Botaurus_minutus" %in% phy$tip.label
phy$tip.label[which(phy$tip.label=="Botaurus_minutus")] = "Ixobrychus_minutus"

# Larus_bulleri = Chroicocephalus_bulleri
"Chroicocephalus_bulleri" %in% phy$tip.label
phy$tip.label[which(phy$tip.label=="Chroicocephalus_bulleri")] = "Larus_bulleri"

# Larus_novaehollandiae =  	Chroicocephalus_novaehollandiae
"Chroicocephalus_novaehollandiae" %in% phy$tip.label
phy$tip.label[which(phy$tip.label=="Chroicocephalus_novaehollandiae")] = "Larus_novaehollandiae"

# Larus_ridibundus =  	Chroicocephalus_ridibundus
"Chroicocephalus_ridibundus" %in% phy$tip.label
phy$tip.label[which(phy$tip.label=="Chroicocephalus_ridibundus")] = "Larus_ridibundus"

# Leiopicus_medius = Dendrocoptes_medius
"Dendrocoptes_medius" %in% phy$tip.label
phy$tip.label[which(phy$tip.label=="Dendrocoptes_medius")] = "Leiopicus_medius"


tree_sub <- keep.tip(phy, spp)
# Calculate distances 
D <- cophenetic.phylo(tree_sub)

# Convert names -> integer matrix positions
i <- match(sp1, rownames(D))
j <- match(sp2, colnames(D))

# Vector of distances in original pair order
phylo_dist <- data.table(pd = D[cbind(i, j)],
                         species_pair=pairs)
res = merge(res, phylo_dist, by='species_pair', all.x=T)




## get FD ------------------------------------------------------------------

avonet = fread("data/AVONET1_BirdLife.csv")
avonet$species = gsub(" ", "_", avonet$Species1)
avosub = avonet[species %in% spp, ]

# select traits to use
input <- avosub[,.(species, Beak.Length_Culmen, Beak.Width, Beak.Depth, Tarsus.Length, Wing.Length, 
                   Kipps.Distance, Secondary1, `Hand-Wing.Index`, Tail.Length, Mass)]
gow <- cluster::daisy(input[, !"species"], metric="gower")
attr(gow, "Labels") <- avosub$species
gow_mat <- as.matrix(gow) # rows and columns have species names!

gower_pairs <- data.table(species_pair = pairs,
                          gower = gow_mat[cbind(match(pair_dt[[1]], rownames(gow_mat)), match(pair_dt[[2]], colnames(gow_mat)))]
                          )

res = merge(res, gower_pairs, by='species_pair', all.x=T)
rm(gow_mat, gower_pairs)





## Attach traits -------------------------------------------------------

res[, sp1:=gsub("\\|.*", "", species_pair), ]
res[, sp2:=gsub(".*\\|", "", species_pair), ]


res[, `:=`(
  habitat_sp1 = avonet$Habitat[match(sp1, avonet$species)], # find sp1 in avonet$species
  habitat_sp2 = avonet$Habitat[match(sp2, avonet$species)]
)]


res[, habitat_same:=as.numeric(habitat_sp1==habitat_sp2), ]
res$sp1 = factor(res$sp1)
res$sp2 = factor(res$sp2)




## Attach occupancy --------------------------------------------------------

# occupancy t1, occupancy change, and ratio
files = dir("data", pattern="processed_occupancy_for_CFs_2_", full.names = T)
dat = lapply(files, readRDS)
dat = rbindlist(dat)
occ = dat[, .(total_occ=sum(mean.psi)), .(datasetID, scientificName, endYear)]
occ[datasetID==5, endYear:=ifelse(endYear==3, 2, 1),]
occ[, occ_change:=total_occ[endYear==2] - total_occ[endYear==1], .(datasetID, scientificName)]
occt1 = occ[endYear==1]

res[, `:=`(
  occ_sp1_t1 = occt1$total_occ[match(paste(dataset_id, sp1), paste(occt1$datasetID, occt1$scientificName))],
  occ_sp2_t1 = occt1$total_occ[match(paste(dataset_id, sp2), paste(occt1$datasetID, occt1$scientificName))],
  occ_change_sp1 = occ$occ_change[match(paste(dataset_id, sp1), paste(occ$datasetID, occ$scientificName))],
  occ_change_sp2 = occ$occ_change[match(paste(dataset_id, sp2), paste(occ$datasetID, occ$scientificName))]
)]

res[, occ_ratio:=(min(c(occ_sp1_t1, occ_sp2_t1))/(max(c(occ_sp1_t1, occ_sp2_t1)))), 
    .(species_pair, dataset_id)]
res[, min_occ:=min(c(occ_sp1_t1, occ_sp2_t1)), 
    .(species_pair, dataset_id)]










## Multimembership models ----------------------------------------------------

res$atlas = factor(res$dataset_id, levels = c(5,6,17,26),
                            labels=c("Czechia", "New York", "New Zealand", "Europe"))
for(a in c(5,6,17,26)){
  mod_dat = res[dataset_id==a,]
  mod_dat = na.omit(mod_dat)
  mod_dat = droplevels(mod_dat)
  mod_dat = unique(mod_dat[, .(delta_cor_median, pd, gower, habitat_same, sp1, sp2)])
  
  pair_members <- paste(as.character(mod_dat$sp1), as.character(mod_dat$sp2), sep = ",")
  
  # Construct sparse membership matrix
  W_species <- lmerMultiMember::weights_from_vector(pair_members)
  dim(W_species)
  
  ### base model -----
  m_base <- lmerMultiMember::lmer(data = mod_dat,
                                delta_cor_median ~ pd + gower + habitat_same +(1 | species),
                                memberships = list(species = W_species),
                                REML=F
  )
  # scaling
  vars <- c("pd","gower")
  mod_dat[, (paste0(vars, "_sc")) :=lapply(.SD, function(x) c(scale(x))), .SDcols = vars]
  m_base_sc <- lmerMultiMember::lmer(data = mod_dat, 
                                     delta_cor_median ~ pd_sc + gower_sc + habitat_same +(1|species), 
                                     memberships = list(species = W_species), 
                                     REML = FALSE)
  
  
  
  ### adding habitat to base model -------------
  
  species_trait = avonet[species %in% c(mod_dat$sp1, mod_dat$sp2), .(species, Habitat)]
  table(species_trait$Habitat)
  
  # pick baseline habitat
  species_trait$Habitat = gsub(" ", "_", species_trait$Habitat)
  species_trait$Habitat = as.factor(species_trait$Habitat)
  species_trait$Habitat <- relevel(species_trait$Habitat, ref = "Forest")
  
  dummy <- model.matrix(~ Habitat, data = species_trait)[, -1, drop = FALSE]
  colnames(dummy) <- gsub("^Habitat", "hab_", colnames(dummy))
  rownames(dummy) <- species_trait$species
  
  pair_dummy_sum <- dummy[mod_dat$sp1, ] + dummy[mod_dat$sp2, ]
  mod_dat <- cbind(mod_dat, as.data.frame(pair_dummy_sum))
  
  pair_members <- paste(as.character(mod_dat$sp1), as.character(mod_dat$sp2), 
                        sep = ",")
  # Construct sparse membership matrix
  W_species <- lmerMultiMember::weights_from_vector(pair_members)
  dim(W_species)
  
  if(a==5){
    m_hab <- lmerMultiMember::lmer(
      delta_cor_median ~ pd + gower + habitat_same +
        hab_Grassland+ hab_Human_Modified+ hab_Marine+ hab_Riverine+ hab_Rock+ 
        hab_Shrubland+ hab_Wetland+ hab_Woodland+
        (1 | species),
      data = mod_dat,
      memberships = list(species = W_species), REML=F
    )
  }
  if(a==6){
    m_hab <- lmerMultiMember::lmer(
      delta_cor_median ~ pd + gower + habitat_same +
        hab_Coastal +hab_Grassland +hab_Human_Modified +hab_Marine +hab_Riverine 
      +hab_Shrubland +hab_Wetland +hab_Woodland +
        (1 | species),
      data = mod_dat,
      memberships = list(species = W_species), REML=F
    )
  }
  if(a==17){
    m_hab <- lmerMultiMember::lmer(
      delta_cor_median ~ pd + gower + habitat_same +
        hab_Coastal+ hab_Grassland+ hab_Human_Modified+ hab_Marine+ hab_Riverine+
        hab_Shrubland+ hab_Wetland+ hab_Woodland+
        (1 | species),
      data = mod_dat,
      memberships = list(species = W_species), REML=F
    )
  }
  if(a==26){
    m_hab <- lmerMultiMember::lmer(
      delta_cor_median ~ pd + gower + habitat_same +
        hab_Coastal+ hab_Grassland+ hab_Human_Modified+ hab_Marine+ hab_Riverine+
        hab_Shrubland+ hab_Wetland+ hab_Woodland+ hab_Rock+
        (1 | species),
      data = mod_dat,
      memberships = list(species = W_species), REML=F
    )
  }
  anova(m_base, m_hab)
  
  
  
  
  ### adding primary lifestyle -----
  species_trait = avonet[species %in% c(mod_dat$sp1, mod_dat$sp2), .(species, Primary.Lifestyle)]
  table(species_trait$Primary.Lifestyle)
  
  # pick baseline habitat
  species_trait$Primary.Lifestyle = gsub(" ", "_", species_trait$Primary.Lifestyle)
  species_trait$Primary.Lifestyle = as.factor(species_trait$Primary.Lifestyle)
  species_trait$Primary.Lifestyle <- relevel(species_trait$Primary.Lifestyle, ref = "Generalist")
  
  dummy <- model.matrix(~ Primary.Lifestyle, data = species_trait)[, -1, drop = FALSE]
  colnames(dummy) <- gsub("^Primary.Lifestyle", "pl_", colnames(dummy))
  rownames(dummy) <- species_trait$species
  
  pair_dummy_sum <- dummy[mod_dat$sp1, ] + dummy[mod_dat$sp2, ]
  mod_dat <- cbind(mod_dat, as.data.frame(pair_dummy_sum))
  
  # Construct sparse membership matrix
  pair_members <- paste(as.character(mod_dat$sp1), as.character(mod_dat$sp2), sep = ",")
  W_species <- lmerMultiMember::weights_from_vector(pair_members)
  dim(W_species)
  
  if(a %in% c(5, 6, 17, 26)){
    m_pl <- lmerMultiMember::lmer(
      delta_cor_median ~ pd + gower + habitat_same +
        pl_Aerial + pl_Aquatic + pl_Insessorial + pl_Terrestrial+
        (1 | species),
      data = mod_dat,
      memberships = list(species = W_species), REML = F
    )
  }
  anova(m_base, m_pl)
  
  
  
  ### adding trophic niche -----
  species_trait = avonet[species %in% c(mod_dat$sp1, mod_dat$sp2), .(species, Trophic.Niche)]
  table(species_trait$Trophic.Niche)
  
  
  # pick baseline habitat
  species_trait$Trophic.Niche = gsub(" ", "_", species_trait$Trophic.Niche)
  species_trait$Trophic.Niche = as.factor(species_trait$Trophic.Niche)
  species_trait$Trophic.Niche <- relevel(species_trait$Trophic.Niche, ref = "Omnivore")
  
  dummy <- model.matrix(~ Trophic.Niche, data = species_trait)[, -1, drop = FALSE]
  colnames(dummy) <- gsub("^Trophic.Niche", "tn_", colnames(dummy))
  rownames(dummy) <- species_trait$species
  
  pair_dummy_sum <- dummy[mod_dat$sp1, ] + dummy[mod_dat$sp2, ]
  mod_dat <- cbind(mod_dat, as.data.frame(pair_dummy_sum))
  
  # Construct sparse membership matrix
  pair_members <- paste(as.character(mod_dat$sp1), as.character(mod_dat$sp2), sep = ",")
  W_species <- lmerMultiMember::weights_from_vector(pair_members)
  dim(W_species)
  
  if(a==5){
    m_tn <- lmerMultiMember::lmer(
      delta_cor_median ~ pd + gower + habitat_same +
        tn_Aquatic_predator + tn_Granivore + tn_Herbivore_aquatic + tn_Invertivore + tn_Vertivore +
        (1 | species),
      data = mod_dat,
      memberships = list(species = W_species), REML=F
    )
  }
  if(a==6){
    m_tn <- lmerMultiMember::lmer(
      delta_cor_median ~ pd + gower + habitat_same +
        tn_Aquatic_predator+ tn_Granivore+ tn_Herbivore_aquatic+ tn_Herbivore_terrestrial+ 
        tn_Invertivore+ tn_Scavenger+ tn_Vertivore+
        (1 | species),
      data = mod_dat,
      memberships = list(species = W_species), REML=F
    )
  }
  if(a==17){
    m_tn <- lmerMultiMember::lmer(
      delta_cor_median ~ pd + gower + habitat_same +
        tn_Aquatic_predator+ tn_Frugivore+ tn_Granivore+ tn_Herbivore_aquatic+ 
        tn_Invertivore+ tn_Nectarivore+ tn_Vertivore+
        (1 | species),
      data = mod_dat,
      memberships = list(species = W_species), REML=F
    )
  }
  if(a==26){
    m_tn <- lmerMultiMember::lmer(
      delta_cor_median ~ pd + gower + habitat_same +
        tn_Aquatic_predator + tn_Granivore + tn_Herbivore_aquatic + tn_Herbivore_terrestrial +
        tn_Invertivore + tn_Vertivore +
        (1 | species),
      data = mod_dat,
      memberships = list(species = W_species), REML=F
    )
  }
  anova(m_base, m_tn)

  
  
  
  ### adding occupancy traits -------------
  mod_dat = res[dataset_id==a,]
  mod_dat = na.omit(mod_dat)
  mod_dat = droplevels(mod_dat)
  # add variables
  mod_dat[, `:=`(
    occ_change_mean = (occ_change_sp1 + occ_change_sp2) / 2   # average trend of the pair
  )] 
  
  mod_dat = unique(mod_dat[, .(delta_cor_median, pd, gower, habitat_same, sp1, sp2, 
                               occ_ratio, occ_sp1_t1, occ_sp2_t1, min_occ,
                               occ_change_mean)])
  vars <- c("pd","gower", "min_occ", "occ_sp1_t1", "occ_sp2_t1", "occ_change_mean")
  mod_dat[, (paste0(vars, "_sc")) :=lapply(.SD, function(x) c(scale(x))), .SDcols = vars]

  # Construct sparse membership matrix
  pair_members <- paste(as.character(mod_dat$sp1), as.character(mod_dat$sp2), sep = ",")
  W_species <- lmerMultiMember::weights_from_vector(pair_members)
  
  m_occ <- lmerMultiMember::lmer(data = mod_dat,
                                 delta_cor_median ~ pd_sc + gower_sc + habitat_same + 
                                  min_occ_sc + occ_ratio + occ_change_mean_sc + 
                                  (1 | species),
      memberships = list(species = W_species), 
      REML=F
    )
  anova(m_base_sc, m_occ)
  broom.mixed::tidy(m_occ)
  r2(m_base_sc)
  r2(m_occ)
 
  
  
  
  
  ### phylo signal in the model -----
  
  # since pd is in the model, test the phylogentic signal in the residual species effect;
  
   best_model <- m_occ
  sp_effects <- ranef(best_model)$species   
  tree_subsub = keep.tip(tree_sub, rownames(sp_effects))
  
  sp_effects <- sp_effects[tree_subsub$tip.label, , drop = FALSE]  # match order to tree tips
  effect_vec <- setNames(sp_effects[, 1], rownames(sp_effects))
  
  # Blomberg's K - tests against a Brownian motion expectation
  (blombergs_k = phytools::phylosig(tree_subsub, effect_vec, method = "K", test = TRUE, nsim = 1000))
  
  # Pagel's lambda - estimates how much the covariance matches phylogeny (0 = none, 1 = full BM)
  (pagel_l = phylosig(tree_subsub, effect_vec, method = "lambda", test = TRUE))
  
 

  
  
  # SAVE models
  saveRDS(list(phylo=tree_subsub, 
               atlas=a,
               models = list(
                m_base=m_base, 
                m_hab=m_hab, 
                m_pl=m_pl, 
                m_tn=m_tn, 
                m_base_sc=m_base_sc,
                m_occ = m_occ
               ),
               blomberg=blombergs_k,
               pagel=pagel_l),
          file=paste0("output/model_phy_sig_output_", a, ".rds"))
  
} ###### end model loop

