# Author: Melanie Tietje
# Email: tietje@fzp.czu.cz
# GitHub: @Eryops1


# Load libraries ----------------------------------------------------------

rm(list=(ls())) # clear workspace
library(data.table)
library(ggplot2)
library(terra)
library(sf)
library(fields)
library(cowplot)
library(gstat)
library(sp)
library(ggExtra)
library(ggpubr)
library(ggpmisc)
library(ggsci)
library(scico)
library(scattermore)
library(ggpointdensity)
library(circlize)

theme_set(theme_classic(base_size = 9)+
            theme(strip.background = element_blank())+
            theme(axis.title = ggtext::element_markdown())+
            theme(axis.title.x = ggtext::element_markdown())+
            theme(axis.title.y = ggtext::element_markdown())+
            theme(legend.title = ggtext::element_markdown()))





# --------------------------------------------------------------------------#
# FIGURE 1 ------------------------------------------------------------------
# --------------------------------------------------------------------------#




## maps --------------------------------------------------------------------

# gonna make those smaller, distribution does not tell much. just use grids.

grid = lapply(paste0("data/processed_occupancy_for_CFs_2_",c(5,6,17,26),".rds"), readRDS)
grid = rbindlist(grid)
occ_sum = grid[, .(mean.psi_sum=sum(mean.psi)), .(datasetID, siteID)]

# reread to single geometry
grid = lapply(paste0("data/all_scales_atlas_", c(5,6,17,26), ".gpkg"), vect)
grid[[2]] = grid[[2]][grid[[2]]$scalingID==2,]
grid[[2]]$scalingID=1 # adjustfor occupancy scales
grid = vect(grid)
grid = grid[grid$scalingID==1,]
grid = merge(grid, occ_sum, by=c("datasetID", "siteID"))

shape_names = paste0("data/", c("europe", "czechia", "nz", "new_york_state"), ".gpkg")
shapes = lapply(shape_names, function(x){st_as_sf(vect(x))})
names(shapes) = c("eu", "cz", "nz", "ny")
lapply(shapes, crs, proj=T)


cex = 8
lw = 0.3

# petrs_fav_proj = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
# ended up not using this


p_grid <- st_as_sf(grid)
(maps = plot_grid(
  plot_grid(
    ggplot(shapes$cz)+
      geom_sf(data=p_grid[p_grid$datasetID==5,], aes(fill=mean.psi_sum), show.legend = F, col=NA)+
      geom_sf(col="black", fill="NA", lwd=lw)+
      theme_void()+
      ggtitle("Czechia")+
      scale_fill_gradient("Species richness", low = 'white', high = '#c90f04')+
      theme(title = element_text(size = cex)),
    ggplot(shapes$ny)+
      geom_sf(data=p_grid[p_grid$datasetID==6,], aes(fill=mean.psi_sum), show.legend = F, col=NA)+
      geom_sf(col="black", fill="NA", lwd=lw)+
      theme_void()+
      ggtitle("New York State")+
      scale_fill_gradient("Species richness", low = 'white', high = '#85bb07')+
      theme(title = element_text(size = cex)),
    ncol=1),
  
  ggplot(shapes$nz)+
    geom_sf(data=p_grid[p_grid$datasetID==17,], aes(fill=mean.psi_sum), show.legend = F, col=NA)+
    geom_sf(col="black", fill="NA", lwd=lw)+
    theme_void()+
    ggtitle("New Zealand")+
    scale_fill_gradient("Species richness", low = 'white', high = '#fecd04')+
    theme(title = element_text(size = cex)),
  
  ggplot(shapes$eu)+
    geom_sf(data=p_grid[p_grid$datasetID==26,], aes(fill=mean.psi_sum), show.legend = F, col=NA)+
    geom_sf(col="black", fill="NA", lwd=lw)+
    theme_void()+
    ggtitle("Europe")+
    scale_fill_gradient("Species richness", low = 'white', high = '#5f89d8')+
    theme(title = element_text(size = cex))
  
  , ncol=3, rel_widths = c(1,1.1,1.2))
)

ggsave(plot = maps, filename = "figures2/maps2.svg", height = 5, width = 5)





# lil elevation map
ele <- geodata::elevation_30s("Czechia", path = "data/environment/elevation/elevation_cz")
ele =mask(ele, p_grid)
plot(ele, col=map.pal("grey", 5))
plot(shapes$cz[1], add=T, col=NA, lwd=2)







## Fig 1 timeline ---------------------------------------------------------------
# number of species and pairs in which years
tmp = readRDS("data/final_analysis_input2.rds")
res = tmp[[1]]
pairwise_res = tmp[[2]]
rm(tmp)

pairs_n = pairwise_res[, .N, .(dataset_id)]
sum(pairs_n$N) 

pairwise_res[, sp1:=gsub("\\|.*", "", species_pair), ]
pairwise_res[, sp2:=gsub(".*\\|", "", species_pair), ]
cz_n = length(unique(c(pairwise_res$sp1[pairwise_res$dataset_id==5], pairwise_res$sp2[pairwise_res$dataset_id==5])))
ny_n = length(unique(c(pairwise_res$sp1[pairwise_res$dataset_id==6], pairwise_res$sp2[pairwise_res$dataset_id==6])))
nz_n = length(unique(c(pairwise_res$sp1[pairwise_res$dataset_id==17], pairwise_res$sp2[pairwise_res$dataset_id==17])))
eu_n = length(unique(c(pairwise_res$sp1[pairwise_res$dataset_id==26], pairwise_res$sp2[pairwise_res$dataset_id==26])))



# Average time span covered: ceiling(mean(years(atlas2))) -
# ceiling(mean(years(atlas1)))
mean(ceiling(mean(2014:2017)) - ceiling(mean(1985:1989)), # CZ
     ceiling(mean(2000:2005)) - ceiling(mean(1980:1985)), # NY
     ceiling(mean(1999:2004)) - ceiling(mean(1969:1979)), # NZ
     ceiling(mean(2013:2017)) - ceiling(mean(1972:1995)) # EU
)

c(ceiling(mean(2014:2017)) - ceiling(mean(1985:1989)), # CZ
  ceiling(mean(2000:2005)) - ceiling(mean(1980:1985)), # NY
  ceiling(mean(1999:2004)) - ceiling(mean(1969:1979)), # NZ
  ceiling(mean(2013:2017)) - ceiling(mean(1972:1995)) # EU
)



tmp = data.table(atlas=rep(c("Czechia", "New York", "New Zealand", "Europe"), each=2),
                 lower = c(1985, 2014, 1980, 2000, 1969, 1999, 1972, 2013),
                 upper = c(1989, 2017, 1985, 2005, 1979, 2004, 1995, 2017))

# add landuse
tmp = rbind(tmp, data.table(atlas="Landuse", 
                            upper = c(1987, 2016, 1983, 2003, 1974, 2002, 1984, 2015), 
                            lower = c(1987, 2016, 1983, 2003, 1974, 2002, 1984, 2015))
)

# add clim data times
files = fread("data/clim_file_names.csv")
tmp = rbind(tmp, data.table(atlas="Climate", 
                            upper = c(1995, 2005, 2017), 
                            lower = c(1979, 1999, 2013))
)


tmp$atlas = factor(tmp$atlas, levels=c("Czechia", "Europe", "New York", "New Zealand", "Landuse", "Climate"))
tmp$atlas2 = factor(tmp$atlas, levels=c("Czechia", "Europe", "New York", "New Zealand", "Landuse", "Climate"),
                    labels = c("biological", "biological", "biological", "biological", "environmental", "environmental"))
s = 4
(timeline = ggplot()+
    geom_linerange(data=tmp[!atlas %in% c("Landuse", "Climate"),], aes(x=atlas, ymin=lower, ymax=upper, col=atlas), size=s, show.legend = F)+
    geom_point(data=tmp[atlas %in% c("Landuse"),], aes(x=atlas, y=lower), col=c("darkgreen"), pch=15, show.legend = F, size=s)+
    geom_linerange(data=tmp[atlas %in% c("Climate"),], aes(x=atlas, ymin=lower, ymax=upper), col=c("purple"), show.legend = F, size=s)+
    scale_color_startrek(guide='none')+
    facet_grid(~atlas2, scales="free", space = "free")+
    #scale_size_continuous(guide='none')+
    labs(x="", y="Years")+
    theme(panel.background = element_blank(),
          axis.text.x = element_text(angle=45, hjust = 1)))


ggsave(paste0("figures2/data_timeline.svg"), height=2.2, width=2, bg = NULL)




# toy phylogeny
library(ape)
tree <- rtree(10)  # Create a random tree with 10 tips
ultrametric_tree <- chronos(tree)  # Make it ultrametric
svg("figures2/toy_phylo.svg", width=3,height=3)
ape::plot.phylo(ultrametric_tree, show.tip.label = FALSE)
dev.off()









## d - counterfactuals -----------------------------------------------------

# species map example
dat = readRDS("data/processed_occupancy_for_CFs_2_5.rds")
sp = unique(dat$scientificName)
sub = dat[scientificName==sp[1] & endYear==3,]

## map occupancy values
grid = vect(paste0("data/all_scales_atlas_5.gpkg"))
grid = grid[grid$scalingID==1,]
grid = merge(grid, sub, by=c("datasetID", "scalingID", "siteID"), all.y=TRUE)

p1 = st_as_sf(grid)
ggplot(p1, aes(fill=mean.psi))+
  geom_sf(col=NA)+
  scale_fill_binned("Occupancy", low="black", high="white")+
  theme_void()+
  # theme(legend.key.width = unit(4,"mm"),
  #       legend.key.height = unit(5,"mm"))
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  geom_sf(data=shapes$cz, fill=NA, col="black", lwd=1)
ggsave("figures2/schem_map.svg", width=3, height=2)













# --------------------------------------------------------------------------#
# FIGURE 5 ------------------------------------------------------------------
# --------------------------------------------------------------------------#

# this is simply a collection of toy examples on how this could work, not real data.

## joint movement ----------------------------------------------------------

a1 = c(1,1,1,1,1,1,1,1,1,0,0,0)
b1 = c(0,0,0,1,1,0,0,1,0,0,0,0)
a2 = c(0,0,0,1,1,1,1,1,1,1,1,1)
b2 = c(0,0,0,0,0,0,1,1,0,0,1,0)

# stats
cor.test(a1,b1,method="s")
cor.test(a2,b2,method="s")

# get jaccard for each grid T1 vs T2 (we replace sites (=rows) here with time)
jacc = data.frame(grid=c(1:12), jaccard=NA)
for(i in 1:length(a1)){
  jacc$jaccard[i] <- vegan::vegdist(rbind(cbind(a1[i], b1[i]), cbind(a2[i], b2[i])), method="jaccard", binary=TRUE)
}
jacc$jaccard = abs(jacc$jaccard-1)
mean(jacc$jaccard)
# note to me: vegan returns DISSIMILARITY. "Jaccard index is computed as
# 2B/(1+B), where B is Bray–Curtis dissimilarity"


# plot 
x = c(1,1,1,2,2,2,3,3,3,4,4,4)
y = c(3,2,1,3,2,1,3,2,1,3,2,1)

par(mfrow=c(1,2), mar=c(1,1,1,1))
plot(x-0.25, y, col=c("white", "brown")[a1+1], pch=20, cex=5, xlim=c(0.5,4.5), ylim=c(0.5,3.5),
     xlab="", ylab="", xaxt="n", yaxt="n",
     main=paste0("rho=", round(cor.test(a1,b1,method="s")$estimate, 2), ", jacc_avg=", round(mean(jacc$jaccard),2)))
abline(h=c(1.5, 2.5), v=c(1.5, 2.5, 3.5))
points(x+0.25, y, col=c("white", "turquoise")[b1+1], pch=20, cex=5)

plot(x-0.25, y, col=c("white", "brown")[a2+1], pch=20, cex=5, xlim=c(0.5,4.5), ylim=c(0.5,3.5),
     xlab="", ylab="", xaxt="n", yaxt="n",
     main=paste0("rho=", round(cor.test(a2,b2,method="s")$estimate, 2)))
abline(h=c(1.5, 2.5), v=c(1.5, 2.5, 3.5))
points(x+0.25, y, col=c("white", "turquoise")[b2+1], pch=20, cex=5)
text(x+0.25,y-0.25,jacc$jaccard[])







## expansion (or contraction) -----------------------------------------------

a1 = c(1,1,1,1,1,1,1,1,0,0,0,0)
b1 = c(0,0,0,1,0,0,0,1,0,0,0,0)
a2 = c(1,1,1,1,1,1,1,1,1,0,1,0)
b2 = c(0,0,0,1,1,0,1,1,0,0,0,0)

# stats
cor.test(a1,b1,method="s")
cor.test(a2,b2,method="s")

# get jaccard for each grid T1 vs T2 (we replace sites (=rows) here with time)
jacc = data.frame(grid=c(1:12), jaccard=NA)
for(i in 1:length(a1)){
  jacc$jaccard[i] <- vegan::vegdist(rbind(cbind(a1[i], b1[i]), cbind(a2[i], b2[i])), method="jaccard", binary=TRUE)
}
jacc$jaccard = abs(jacc$jaccard-1)
jacc$jaccard[is.na(jacc$jaccard)] = 1 # replace the no measure with 1 because nothing is changing...
mean(jacc$jaccard, na.rm=TRUE)
# note to me: vegan returns DISSIMILARITY. "Jaccard index is computed as
# 2B/(1+B), where B is Bray–Curtis dissimilarity"


# plot 
par(mfrow=c(1,2), mar=c(1,1,1,1))
plot(x-0.25, y, col=c("white", "brown")[a1+1], pch=20, cex=5, xlim=c(0.5,4.5), ylim=c(0.5,3.5),
     xlab="", ylab="", xaxt="n", yaxt="n",
     main=paste0("rho=", round(cor.test(a1,b1,method="s")$estimate, 2), ", jacc_avg=", round(mean(jacc$jaccard, na.rm=T),2)))
abline(h=c(1.5, 2.5), v=c(1.5, 2.5, 3.5))
points(x+0.25, y, col=c("white", "turquoise")[b1+1], pch=20, cex=5)

plot(x-0.25, y, col=c("white", "brown")[a2+1], pch=20, cex=5, xlim=c(0.5,4.5), ylim=c(0.5,3.5),
     xlab="", ylab="", xaxt="n", yaxt="n",
     main=paste0("rho=", round(cor.test(a2,b2,method="s")$estimate, 2)))
abline(h=c(1.5, 2.5), v=c(1.5, 2.5, 3.5))
points(x+0.25, y, col=c("white", "turquoise")[b2+1], pch=20, cex=5)
text(x+0.25,y-0.25,jacc$jaccard[])








## occupancy change --------------------------------------------------------

a1 = c(2,2,1,3,3,2,2,2,1,0,0,0)
a2 = c(2,0,1,3,0,2,2,2,1,0,0,0)
b1 = c(0,0,0,1,1,0,1,1,0,0,0,0)
b2 = c(0,0,0,2,0,0,1,1,0,0,0,0)

# stats
cor.test(a1,b1,method="s")
cor.test(a2,b2,method="s")

# get jaccard for each grid T1 vs T2 (we replace sites (=rows) here with time)
jacc = data.frame(grid=c(1:12), jaccard=NA)
for(i in 1:length(a1)){
  jacc$jaccard[i] <- vegan::vegdist(rbind(cbind(a1[i], b1[i]), cbind(a2[i], b2[i])), method="jaccard", binary=TRUE)
}
jacc$jaccard = abs(jacc$jaccard-1)
jacc$jaccard[is.na(jacc$jaccard)] = 1 # replace the no measure with 1, because nothing is changing...
mean(jacc$jaccard, na.rm=TRUE)
# note to me: vegan returns DISSIMILARITY. "Jaccard index is computed as
# 2B/(1+B), where B is Bray–Curtis dissimilarity"


# plot 
par(mfrow=c(1,2), mar=c(1,1,1,1))
plot(x-0.25, y, col="brown", pch=20, cex=c(0,2,4,6)[a1+1], xlim=c(0.5,4.5), ylim=c(0.5,3.5),
     xlab="", ylab="", xaxt="n", yaxt="n",
     main=paste0("rho=", round(cor.test(a1,b1,method="s")$estimate, 2), ", jacc_avg=", round(mean(jacc$jaccard),2)))
abline(h=c(1.5, 2.5), v=c(1.5, 2.5, 3.5))
points(x+0.25, y, col="turquoise", pch=20, cex=c(0,2,4,6)[b1+1])

plot(x-0.25, y, col="brown", pch=20, cex=c(0,2,4,6)[a2+1], xlim=c(0.5,4.5), ylim=c(0.5,3.5),
     xlab="", ylab="", xaxt="n", yaxt="n",
     main=paste0("rho=", round(cor.test(a2,b2,method="s")$estimate, 2)))
abline(h=c(1.5, 2.5), v=c(1.5, 2.5, 3.5))
points(x+0.25, y, col="turquoise", pch=20, cex=c(0,2,4,6)[b2+1])
text(x+0.25,y-0.25,jacc$jaccard[])


# if this was simply binary, correlation would change:

a1 = c(0,1)[as.numeric(a1>0)+1]
a2 = c(0,1)[as.numeric(a2>0)+1]
b1 = c(0,1)[as.numeric(b1>0)+1]
b2 = c(0,1)[as.numeric(b2>0)+1]

# stats
cor.test(a1,b1,method="s")
cor.test(a2,b2,method="s")

# get jaccard for each grid T1 vs T2 (we replace sites (=rows) here with time)
jacc = data.frame(grid=c(1:12), jaccard=NA)
for(i in 1:length(a1)){
  jacc$jaccard[i] <- vegan::vegdist(rbind(cbind(a1[i], b1[i]), cbind(a2[i], b2[i])), method="jaccard", binary=TRUE)
}
jacc$jaccard = abs(jacc$jaccard-1)
jacc$jaccard[is.na(jacc$jaccard)] = 1 # replace the no measure with 1, because nothing is changing...
mean(jacc$jaccard, na.rm=TRUE)
# note to me: vegan returns DISSIMILARITY. "Jaccard index is computed as
# 2B/(1+B), where B is Bray–Curtis dissimilarity"


# plot 
par(mfrow=c(1,2), mar=c(1,1,1,1))
plot(x-0.25, y, col="brown", pch=20, cex=c(0,2,4,6)[a1+1], xlim=c(0.5,4.5), ylim=c(0.5,3.5),
     xlab="", ylab="", xaxt="n", yaxt="n",
     main=paste0("rho=", round(cor.test(a1,b1,method="s")$estimate, 2), ", jacc_avg=", round(mean(jacc$jaccard),2)))
abline(h=c(1.5, 2.5), v=c(1.5, 2.5, 3.5))
points(x+0.25, y, col="turquoise", pch=20, cex=c(0,2,4,6)[b1+1])

plot(x-0.25, y, col="brown", pch=20, cex=c(0,2,4,6)[a2+1], xlim=c(0.5,4.5), ylim=c(0.5,3.5),
     xlab="", ylab="", xaxt="n", yaxt="n",
     main=paste0("rho=", round(cor.test(a2,b2,method="s")$estimate, 2)))
abline(h=c(1.5, 2.5), v=c(1.5, 2.5, 3.5))
points(x+0.25, y, col="turquoise", pch=20, cex=c(0,2,4,6)[b2+1])
text(x+0.25,y-0.25,jacc$jaccard[])




# find more solutions
set.seed(3)
a1 = c(2,2,1,3,3,2,2,2,1,0,0,0)
perm <- sample(c(1:12),replace = F)
Bvals <- seq(0,1,length.out=12)
b1 <- Bvals[perm]
rho_target <- cor(a1, b1, method="spearman")
# build a2. decrease the lowest ranks
ord <- order(a1)
zero_idx <- ord[1:5] # zero the 5 lowest-ranked a1 values
a2 <- a1
a2[zero_idx] <- 0
a2[-zero_idx] <- a2[-zero_idx] * 0.3

# permute b2 to match rho target
nperm <- 200000
best_diff <- Inf
best_b2 <- NULL

for(i in 1:nperm){
  b2_try <- sample(b1, length(b1), replace=FALSE)
  r <- cor(a2, b2_try, method="spearman")
  diff <- abs(r - rho_target)
  if(diff < best_diff){
    best_diff <- diff
    best_b2 <- b2_try
    if(diff < 0.005) break
  }
  cat(i, "\r")
}
b2 <- best_b2

cor(a1, b1, method = "spearman")
cor(a2, b2, method = "spearman")

# get jaccard for each grid T1 vs T2 (we replace sites (=rows) here with time)
jacc = data.frame(grid=c(1:12), jaccard=NA)
for(i in 1:length(a1)){
  jacc$jaccard[i] <- vegan::vegdist(rbind(cbind(a1[i], b1[i]), cbind(a2[i], b2[i])), method="jaccard", binary=TRUE)
}
jacc$jaccard = abs(jacc$jaccard-1)
jacc$jaccard[is.na(jacc$jaccard)] = 1 # replace the no measure with 1, because nothing is changing...
mean(jacc$jaccard, na.rm=TRUE)
# note to me: vegan returns DISSIMILARITY. "Jaccard index is computed as
# 2B/(1+B), where B is Bray–Curtis dissimilarity"


# plot 
par(mfrow=c(1,2), mar=c(1,1,1,1))
plot(x-0.25, y, col="brown", pch=20, cex=a1, xlim=c(0.5,4.5), ylim=c(0.5,3.5),
     xlab="", ylab="", xaxt="n", yaxt="n",
     main=paste0("rho=", round(cor.test(a1,b1,method="s")$estimate, 2), ", jacc_avg=", round(mean(jacc$jaccard),2)))
abline(h=c(1.5, 2.5), v=c(1.5, 2.5, 3.5))
points(x+0.25, y, col="turquoise", pch=20, cex=b1)

plot(x-0.25, y, col="brown", pch=20, cex=a2, xlim=c(0.5,4.5), ylim=c(0.5,3.5),
     xlab="", ylab="", xaxt="n", yaxt="n",
     main=paste0("rho=", round(cor.test(a2,b2,method="s")$estimate, 2)))
abline(h=c(1.5, 2.5), v=c(1.5, 2.5, 3.5))
points(x+0.25, y, col="turquoise", pch=20, cex=b2)
text(x+0.25,y-0.25,jacc$jaccard[])


















# -----------------------------------------------------------------------------#
# MAIN FIGURES ----------------------------------------------------------------
# -----------------------------------------------------------------------------#

rm(list = ls())
source("99_functions.R")

analysis_input <- readRDS("data/final_analysis_input2.rds")
res <- analysis_input[[1]]
pairwise_res <- analysis_input[[2]]
rm(analysis_input)



## FIGURE 2-------------------------------------------------------------------
### caterpillar / histograms and circlize + bray curtis -----------------------

pairwise_res$atlas <- factor(pairwise_res$dataset_id)
pairwise_res$atlas = factor(pairwise_res$dataset_id, levels = c(5,6,17,26),
                            labels=c("Czechia", "New York", "New Zealand", "Europe"))
pairwise_res$atlas = factor(pairwise_res$atlas,
                            levels=c("Czechia", "Europe", "New York", "New Zealand"))

res$atlas = factor(res$dataset_id, levels = c(5,6,17,26),
                   labels=c("Czechia", "New York", "New Zealand", "Europe"))
res$atlas = factor(res$atlas,
                   levels=c("Czechia", "Europe", "New York", "New Zealand"))


(p_cat <- ggplot(pairwise_res[treatment=="observed"],
                 aes(y=y_key, col=greater_less=="neutral"))+
    geom_linerange(data=pairwise_res[treatment!="AM" & greater_less=="neutral", ], 
                   aes(xmin=lower025, xmax=upper975), lwd=0.3, alpha=0.1, show.legend = F, col="grey")+
    geom_linerange(data=pairwise_res[treatment!="AM" & greater_less=="neutral", ], 
                   aes(xmin=lower, xmax=upper))+
    geom_linerange(data=pairwise_res[treatment!="AM" & greater_less!="neutral", ], 
                   aes(xmin=lower025, xmax=upper975), lwd=0.3, alpha=0.1, show.legend = F, col="grey")+
    geom_linerange(data=pairwise_res[treatment!="AM" & greater_less!="neutral", ], 
                   aes(xmin=lower, xmax=upper))+
    scale_color_scico_d("95% quantile uncertainty range overlap 0", 
                        palette=scico_palette_names(categorical = FALSE)[25], begin=0.3, end=0.7)+
    facet_wrap(~atlas, scales="free_y", ncol=4)+
    geom_vline(xintercept = 0, col="black", lty=2)+
    labs(y="species pairs", x="\U0394\U03C1")+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = "bottom", legend.direction = "horizontal"))
ggsave(plot=p_cat, "figures2/caterpillar.png", width = 7, height=7, dpi=300)



atlas_avg = pairwise_res[treatment=="observed", .(atlas_median=median(delta_cor_median)), .(dataset_id, atlas)]
atlas_avg2 = res[treatment=="observed", .(atlas_median=median(delta_cor)), .(dataset_id, atlas)]
(p_pop_it2 = ggplot(res[treatment=="observed"])+
    geom_density(aes(x=delta_cor, fill=factor(atlas), group=iteration), alpha=0.01, col=NA, show.legend = F)+
    geom_density(aes(x=delta_cor, fill=NA), 
                 lty=2, show.legend = T, col="grey20")+
    geom_density(data=pairwise_res[treatment=="observed"], aes(x=delta_cor_median, fill=NA), 
                 lty=1, show.legend = T, col="grey20")+
    labs(x="\U0394\U03C1", y="Frequency")+
    facet_wrap(~atlas, scales="free_y", ncol=4)+
    geom_vline(data=atlas_avg2, aes(xintercept=atlas_median, col=factor(atlas)), lty=1, col="grey20")+
    #geom_vline(data=atlas_avg, aes(xintercept=atlas_median, col=factor(atlas)), lty=1, col="grey20")+
    scale_y_sqrt()+
    scale_fill_startrek(guide="none")+
    scale_color_startrek(guide="none")
)

# IQR
pairwise_res[order(atlas), .(q1=quantile(delta_cor_median, c(0.25)),
                             q3=quantile(delta_cor_median, c(0.75))), .(atlas, treatment)]
pairwise_res[, .(q1=quantile(delta_cor_median, c(0.25)),
                 q3=quantile(delta_cor_median, c(0.75))), .(treatment)]





## FIGURE barplot ## != ZERO -------------------------------------------------

# stats table
stats_tab = pairwise_res[treatment!="AM", .N, 
                         .(atlas, treatment, p.greater_zero<0.05, p.less_zero<0.05)][, perc := 100 * N / sum(N),  
                                                                                     .(atlas, treatment)]
stats_tab_pos = stats_tab[p.greater_zero!=T&p.less_zero!=T,] 
knitr::kable(stats_tab_pos[order(treatment, atlas)], format="rst", digits = 2)

stats_tab[order(treatment), larger_less:=ifelse(p.greater_zero==T, "\U0394\U03C1 > 0", 
                                                ifelse(p.less_zero==T, "\U0394\U03C1 < 0", "neutral"))]
stats_tab$larger_less = factor(stats_tab$larger_less, levels = c("\U0394\U03C1 > 0",
                                                                 "neutral",
                                                                 "\U0394\U03C1 < 0"))

stats_tab$treatment = factor(stats_tab$treatment, levels=c("observed", "S1", "S2", "S3", "S4"))
ggplot(stats_tab, aes(y=perc, x=larger_less, fill=treatment, label=round(perc)))+
  geom_bar(stat = "identity", show.legend = F, alpha=1)+
  geom_label(show.legend = F, fill = NA, border.colour = NA, nudge_y=7, size = 2.5)+
  #scale_y_sqrt(limits=c(0,120), labels=c(0,5,25,50,75,100))+
  facet_grid(atlas~treatment)+
  labs(x="", y="% species pairs")+
  scale_fill_scico_d(palette=scico_palette_names(categorical = FALSE)[25], begin = c(0.9), end=0.3)+
  theme(axis.text.x = element_text(angle=45, hjust = 1), strip.text.x.top= element_text(size=6))
ggsave("figures2/diff_from_zero.png", width = 7, height=5, dpi=300)








## FIGURE 3-------------------------------------------------------------------
### observed vs counterfactual scenarios  -------------------------------------

# stats
res = res[treatment!="AM",]
res = droplevels(res)
res_wide <- dcast(res,
                  species_pair + atlas + iteration ~ treatment,
                  value.var = "delta_cor")

compare = as.character(unique(res$treatment))

# get the differences for each iteration and treatment: observed - treatment
diffs <- rbindlist(lapply(compare, function(tr) {
  d <- res_wide[, .(species_pair, atlas, iteration, delta_diff = observed - get(tr))]
  d[, treatment := tr]
  d
}))

# get N iterations where delta diff is smaller or larger zero (i.e. p-value)
diffs[, `:=`(
  p_greater = (sum(delta_diff <= 0, na.rm=T) + 1) / (.N + 1),  
  p_less    = (sum(delta_diff >= 0, na.rm=T) + 1) / (.N + 1)
), .(species_pair, atlas, treatment)]
# p_greater = fraction of iterations where delta diff < 0 (i.e. obs < sim)
# p_less = fraction of iterations where delta diff > 0 (i.e. obs > sim)
#
# p_greater < 0.05 → obs > sim significantly → "sim < obs"
# p_less < 0.05 → obs < sim significantly → "sim > obs"

diffs_sing=unique(diffs[, .(species_pair, atlas, treatment, p_greater, p_less)])
sig_level = 0.05
diffs_sing[ ,larger_less_sim:=ifelse(p_greater<sig_level, "sim < obs", 
                                     ifelse(p_less<sig_level, "sim > obs", "neutral"))]

pairwise_res = merge(pairwise_res, diffs_sing, all.x=T)

tmp = pairwise_res[!treatment %in% c("AM", "observed"),]
cols = scico(2, palette="lipari", begin=0.3, end=0.7, direction = -1)
#median(tmp$delta_cor_obs[dataset_id==5 & status=="sim > obs"])
(p_scat= ggplot(tmp, 
                aes(y=delta_cor_obs, x=delta_cor_median))+
    geom_abline(slope = c(0, 1e+10, 1), lty=2, col="grey")+
    #geom_pointdensity(alpha=0.1)+ 
    geom_point(data=tmp[larger_less_sim=="neutral"], alpha=0.05, aes(col=larger_less_sim))+
    geom_point(data=tmp[larger_less_sim!="neutral"], alpha=0.05, aes(col=larger_less_sim))+
    scale_color_manual("", values = c(cols[2], cols[1], cols[1]))+
    stat_cor(method = "spearman", size=2.5, cor.coef.name = "rho", aes(label=..r.label..),
             col="black", label.y = +1, label.x = -0.92, r.accuracy = 0.01, p.accuracy = 0.001)+
    #stat_ma_line(col="darkgreen")+ # this uses lmodel2
    coord_cartesian(
      #      ylim=c(round(range(tmp$delta_cor_median)[1]),round(range(tmp$delta_cor_median)[2])), 
      ylim=c(-1,1),
      xlim=c(-1,1))+
    #      xlim=c(round(range(tmp$delta_cor_median)[1]),round(range(tmp$delta_cor_median)[2])))+
    labs(y="Observed median \U0394\U03C1", x="Simulated median \U0394\U03C1")+
    facet_grid(atlas~treatment)+
    theme(strip.text.y.right = element_text(angle = -90),
          legend.margin = margin(0,0,0,0,"mm"),
          legend.text = element_text(margin=margin(0,0,0,0,"mm")), 
          legend.box.spacing = margin(0,0,0,0,"mm"))+
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
)

# violin plot
pairwise_res$treatment = factor(pairwise_res$treatment, 
                                levels=c("AM","S4","S3","S2","S1", "observed"))
(p_box = ggplot(pairwise_res[treatment!="AM"], aes(x=delta_cor_median, y=treatment, fill=treatment))+
    geom_rect(aes(ymin = c(1.5), ymax =c(2.5), xmin = -Inf, xmax = Inf), 
              fill="grey96", col=NA)+
    geom_rect(aes(ymin = c(3.5), ymax =c(4.5), xmin = -Inf, xmax = Inf), 
              fill="grey96", col=NA)+
    # geom_rect(aes(ymin = c(5.5), ymax =c(6.5), xmin = -Inf, xmax = Inf), 
    #           fill="grey94", col=NA)+
    geom_violin(show.legend = F, draw_quantiles = c(0.5), lwd=0.3, alpha=0.9)+
    labs(x="Median \U0394\U03C1", y="")+
    scale_fill_scico_d(palette=scico_palette_names(categorical = FALSE)[25], begin = c(0.3))+
    theme(legend.position = "inside", legend.position.inside = c(0.9,0.8),
          axis.ticks.y=element_blank())+
    facet_wrap("atlas", ncol=1)+
    geom_vline(xintercept = 0, lty=2)
)


(p_hist_sim_obs = ggplot(unique(tmp), aes(x=abs(delta_cor_median-delta_cor_obs)))+
    #geom_abline(slope = c(0, 1e+10, 1), lty=2, col="grey")+
    geom_histogram(alpha=1, aes(fill=larger_less_sim!="neutral"))+
    stat_central_tendency(geom="line", lwd=0.4, type="median")+
    scale_fill_manual("", values = c(cols[2], cols[1], cols[1]), labels = c("neutral", "sim\U2260obs"))+
    labs(x="|Simulated - observed median \U0394\U03C1|", y="Species pairs count")+
    facet_grid(atlas~treatment, scales="free_y")+
    theme(strip.text.y.right = element_text(angle = -90),
          legend.margin = unit(c(0,0,0,0),"mm"))+
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
)




plot_grid(p_box, p_scat, labels = "auto", rel_widths = c(1,1.7))
ggsave("figures2/box_scatter.png", width=9, height = 5, dpi=300)







## add cat sim-/obs-, sim+/obs+
# pairwise_res[, story := fcase(
#   delta_cor_obs > 0 & delta_cor_median > 0, "obs + / sim +",
#   delta_cor_obs < 0 & delta_cor_median < 0, "obs - / sim -",
#   delta_cor_obs > 0 & delta_cor_median < 0, "obs + / sim -",
#   delta_cor_obs < 0 & delta_cor_median > 0, "obs - / sim +",
#   default = "near zero / ambiguous"
# )]
# pairwise_res$story = factor(pairwise_res$story, levels = c("obs + / sim +", "obs + / sim -",  "obs - / sim +",  "obs - / sim -"))
# 
# 
# quadrant_stats <- pairwise_res[order(atlas, treatment), 
#                                .N, .(atlas, treatment, story, larger_less_sim)][
#                                  , perc := 100 * N / sum(N, na.rm=T),  .(atlas, treatment)]
# quadrant_stats = quadrant_stats[!treatment %in% c("observed", "AM")]
# quadrant_stats <- droplevels(quadrant_stats)
# 
# (p_bar_groups = ggplot(quadrant_stats,  aes(x=story, y=perc, fill=story, 
#                                             alpha=larger_less_sim!="neutral", label=round(perc,1)))+
#     geom_bar(stat="identity", show.legend = F, col="black", lwd=0.1)+
#     facet_grid(atlas~treatment, scales="free_y")+
#     labs(x="", y="% species pairs")+
#     #scale_fill_scico_d(palette=scico_palette_names(categorical = FALSE)[25], begin = c(0.3), end=0.7)+
#     scale_fill_manual(values=c("darkred", "darkblue", "yellow3","darkgreen"))+
#     theme(axis.text.x = element_text(angle=50, hjust=1))
# )
# 
# 
# pairwise_res[, scat_cat:=ifelse(larger_less_sim=="neutral", "neutral", as.character(story))]
# pairwise_res$scat_cat = factor(pairwise_res$scat_cat, levels = c("neutral", 
#                                                                  "obs + / sim +",
#                                                                  "obs + / sim -",
#                                                                  "obs - / sim +",
#                                                                  "obs - / sim -"))
# tmp = pairwise_res[!treatment%in%c("observed", "AM")]
# (p_scat2= ggplot(tmp, 
#                  aes(y=delta_cor_obs, x=delta_cor_median))+
#     geom_abline(slope = c(0, 1e+10, 1), lty=2, col="grey70")+
#     geom_point(data=tmp[scat_cat=="neutral"], alpha=0.01, aes(col=scat_cat))+
#     geom_point(data=tmp[scat_cat!="neutral"], alpha=0.05, aes(col=scat_cat))+
#     scale_color_manual(values=c("grey", "darkred", "darkblue", "yellow3","darkgreen"))+
#     stat_cor(method = "spearman", size=2.5, cor.coef.name = "rho", aes(label=..r.label..),
#              col="black", label.y = -1, label.x = -0.9, r.accuracy = 0.01, p.accuracy = 0.001)+
#     #stat_ma_line(col="darkgreen")+ # this uses lmodel2
#     coord_cartesian(ylim=c(round(range(tmp$delta_cor_median)[1]),round(range(tmp$delta_cor_median)[2])), 
#                     xlim=c(round(range(tmp$delta_cor_median)[1]),round(range(tmp$delta_cor_median)[2])))+
#     labs(y="Observed median \U0394\U03C1", x="Simulated median \U0394\U03C1")+
#     facet_grid(dataset_id~treatment)+
#     theme(strip.text.y.right = element_text(angle = 0), legend.position = "none")+
#     guides(colour = guide_legend(override.aes = list(alpha = 1)))
# )
# 
# plot_grid(ncol=2, labels = "auto", p_box, p_scat, rel_widths = c(1,1.6))
# ggsave("figures2/CF_comp_v1.png", width=10, height = 5, dpi=150)
# 
# plot_grid(ncol=3, labels = "auto", p_box, p_scat2, p_bar_groups, rel_widths = c(1, 1.6,0.7))
# ggsave("figures2/CF_comp_v2.png", width=10, height = 5, dpi=150)
# 
# 
# 


# show observed change distributions of cases where simulation does over or underpredict!
pairwise_res[, larger_less_sim_median:=median(delta_cor_obs), .(atlas, treatment, larger_less_sim)]

p_which_differ = ggplot(pairwise_res[!treatment%in%c("observed", "AM") & larger_less_sim!="neutral"], 
                        aes(x=delta_cor_obs, fill=larger_less_sim))+ # larger_less_sim
  geom_vline(xintercept=0, lty=2)+
  geom_histogram(show.legend=T, alpha=1)+
  geom_vline(aes(xintercept=larger_less_sim_median, col=larger_less_sim), show.legend = T, lty=1)+
  scale_fill_scico_d("", palette=scico_palette_names(categorical = FALSE)[25], begin = c(0), end=0.9)+
  scale_color_scico_d("", palette=scico_palette_names(categorical = FALSE)[25], begin = c(0), end=0.9)+
  labs(y="Species pairs count", x="Observed median \U0394\U03C1")+
  facet_grid(atlas~treatment, scales="free_y")



plot_grid(p_hist_sim_obs, p_which_differ, ncol = 1, labels = "auto")
ggsave("figures2/p_hist_sim_obs.png", width=8, height = 7, dpi=300)








## Qualitative transitions -----------------------------------------------------------

tmp = readRDS("data/qualitative_transition_analysis.rds")
list2env(tmp, envir = .GlobalEnv)


# atli
atli = sort(unique(circ_in$atlas))
grid.col = c("segregated"="#1398E9", "neutral"="grey", "aggregated"="#E96513")


png("figures2/circlize.png", width =11, height=11.4, 
    units = "cm", bg = "transparent", res = 300)
par(mfrow=c(2,2))
for(a in 1:length(atli)){
  tmp_iter = circ_in[atlas==atli[a]]
  
  tmp = tmp_iter[, .(median = median(perc), q025 = quantile(perc, 0.025), q975 = quantile(perc, 0.975)),
                 by = .(transition, from, to)]
  tmp_p = tmp[,.(from, to, median)]
  
  circos.par(canvas.xlim = c(-1.05, 1.05), canvas.ylim = c(-1.05, 1.1))
  chordDiagram(tmp_p,
               directional = 1, direction.type = c("arrows"),
               link.arr.type = "big.arrow",
               link.arr.length = 0.1,
               self.link = 1,
               link.largest.ontop = FALSE, grid.col = grid.col,
               transparency = 0.2,
               link.target.prop = TRUE,
               annotationTrack = c("grid"),
               annotationTrackHeight = 0.05,
               scale=F,
               order = c("aggregated", "neutral", "segregated"))
  
  # sum WITHIN each iteration first, then summarize ACROSS iterations
  total_iter   = tmp_iter[, .(val = sum(perc)), by = .(from, iteration)]
  leaving_iter = tmp_iter[from != to, .(val = sum(perc)), by = .(from, iteration)]
  
  summarize_stat = function(dt) {
    dt[, .(median = median(val), q025 = quantile(val, 0.025), q975 = quantile(val, 0.975)), by = from]
  }
  
  total_stats   = summarize_stat(total_iter)
  leaving_stats = summarize_stat(leaving_iter)
  setnames(total_stats,   c("median","q025","q975"), paste0("total_", c("median","q025","q975")))
  setnames(leaving_stats, c("median","q025","q975"), paste0("leaving_", c("median","q025","q975")))
  stats_dt = merge(total_stats, leaving_stats, by = "from", all.x = TRUE)
  stats_dt[is.na(stats_dt)] = 0
  
  cat_stats = split(stats_dt[, -"from"], stats_dt$from)
  cat_stats = lapply(cat_stats, as.list)
  
  fmt = function(m, lo, hi) paste0(format(round(m,1), nsmall=1), "% [", format(round(lo,1), nsmall=1), "-", format(round(hi,1), nsmall=1), "]")
  
  circos.track(ylim = c(0, 1), track.height = 0.25,
               panel.fun = function(x, y) {
                 sector = get.cell.meta.data("sector.index")
                 x_pos = mean(get.cell.meta.data("xlim"))
                 
                 if (sector %in% names(cat_stats)) {
                   s = cat_stats[[sector]]
                   
                   circos.text(x = x_pos, y = 2,
                               labels = paste0(sector, ": ", fmt(s$total_median, s$total_q025, s$total_q975)),
                               facing = "bending.inside", niceFacing = TRUE,
                               cex = 0.7, col = "black", font=1)
                   circos.text(x = x_pos, y = 1.6,
                               labels = paste0("changing: ", fmt(s$leaving_median, s$leaving_q025, s$leaving_q975)),
                               facing = "bending.inside", niceFacing = TRUE,
                               cex = 0.55, col = "grey30", font=3)}},
               bg.border = NA)
  title(main = atli[a], cex.main = 1, font.main = 1)
  circos.clear()
}

dev.off()










## Per species Bray Curtis dissimilarity ---------------------------------
# Step 1: get BrayCurtis matrices (this is similar to Spearman correlation, but
# avoids negative values)

#res = readRDS("data/trait_analysis_input2.rds")
#res = res[treatment=="uncertainty"]

atlas_list = list()
atli = c(5,6,17,26)
mat = list()
for(a in 1:length(atli)){
  dat = readRDS(paste0("data/processed_occupancy_for_CFs_2_", atli[a], ".rds"))
  #load also all the simulation matrices
  species = unique(dat$scientificName)
  files = dir("data/uncertainty", full.names = T, 
              pattern = paste0(paste0(paste0("simulated.*", species, "_", atli[a]), collapse = "|")))
  sims = lapply(files, readRDS)
  names(sims) = species
  sims_t1 = sapply(sims, "[", 1)
  sims_t1 = lapply(sims_t1, as.data.table)
  sims_t2 = sapply(sims, "[", max(lengths(sims)))
  sims_t2 = lapply(sims_t2, as.data.table)
  replacement_t1 = rbindlist(sims_t1, idcol = "scientificName")
  replacement_t2 = rbindlist(sims_t2, idcol = "scientificName")
  
  times = sort(unique(dat$endYear))
  mat_int_t1 = list()
  mat_int_t2 = list()
  for(s in 1:100){
    # processing t1
    tmp <- replacement_t1[, .(scientificName, get(paste0("V", s)))]
    tmp$siteID = rep(c(1:nrow(sims[[1]][[1]])), length(unique(tmp$scientificName)))
    tmp2 <-dcast(tmp, siteID ~ scientificName, value.var = "V2") # this is always V2, thats fine
    tmp2 = vegan::vegdist(t(tmp2[, -c("siteID")]), method="bray", upper=TRUE, diag=TRUE)
    mat_int_t1[[s]] = list(data.table(dataset_id = atli[a],
                                      time_bin = 1), 
                           tmp2)
    # processing t2
    tmp <- replacement_t2[, .(scientificName, get(paste0("V", s)))]
    tmp$siteID = rep(c(1:nrow(sims[[1]][[1]])), length(unique(tmp$scientificName)))
    tmp2 <-dcast(tmp, siteID ~ scientificName, value.var = "V2") # this is always V2, thats fine
    tmp2 = vegan::vegdist(t(tmp2[, -c("siteID")]), method="bray", upper=TRUE, diag=TRUE)
    mat_int_t2[[s]] = list(data.table(dataset_id = atli[a],
                                      time_bin = 2), 
                           tmp2)
  }
  mat[[a]] = list(mat_int_t1, mat_int_t2)
}
names(mat) = atli

# Step 2: get BrayCurtis comparing BrayCurtis dissimilarity with other species
# across sampling periods, for each focal species. Think about it as comparing
# the composition of spatial overlaps that a species has - the group of friends it has!

bc_res = list()
for(a in 1:length(atli)){
  time1 = mat[[a]][[1]] # 100 lists with [1] dt and [2] matrix, T1
  time2 = mat[[a]][[2]] # 100 lists with [1] dt and [2] matrix, T2
  sim_list <- list()
  for(s in 1:100){
    time1_s = time1[[s]][[2]]
    time2_s = time2[[s]][[2]]
    tmp = data.table(species=colnames(time1_s),
                     bray_curtis = NA, 
                     atlas=atli[a])
    time1_m = as.matrix(time1_s)
    time2_m = as.matrix(time2_s)
    for(l in 1:ncol(time1_m)){
      bc = vegan::vegdist(t(data.frame(time1_m[,l], time2_m[,l])))
      tmp$bray_curtis[l] = as.numeric(bc)
    }
    tmp$iteration=s
    sim_list[[s]] = tmp
    
  }
  sim_dt = rbindlist(sim_list)
  bc_res[[a]] = sim_dt
}
rm(sim_list, sims, sims_t1, sims_t2, time1, time2, time1_m, time2_m)
gc()

# now we got 100 values per species personal group

bc = rbindlist(bc_res)
#saveRDS(bc, "data/braycurtis_sims.rds")
#bc = readRDS("data/braycurtis_sims.rds")
bc$atlas = factor(bc$atlas, levels = c(5,6,17,26),
                  labels=c("Czechia", "New York", "New Zealand", "Europe"))
bc$atlas = factor(bc$atlas, levels=c("Czechia", "Europe", "New York", "New Zealand"))

bc_avg = bc[, .(bc_median = median(bray_curtis),
                bc_mad = mad(bray_curtis),
                q025=quantile(bray_curtis, 0.025),
                q975=quantile(bray_curtis, 0.975)), .(atlas, species)]
bc_avg2 <- bc_avg[, .(atlas_avg=median(bc_median)), .(atlas)]


(bcfig = ggplot(bc)+
    geom_density(aes(x=bray_curtis, fill=factor(atlas), group=species), alpha=0.05, col=NA, show.legend = F)+
    geom_density(data=bc_avg, aes(x=bc_median, fill=NA), alpha=1, show.legend = T, col="grey20")+
    geom_vline(data=bc_avg2, aes(xintercept=atlas_avg), lty=1, col="grey20")+
    labs(x="Temporal Bray-Curtis dissimilarity per species", y="Frequency")+
    facet_wrap(~atlas, scales="free_y")+
    scale_y_sqrt()+
    scale_fill_startrek(guide="none")+
    scale_color_startrek(guide="none")
)

ggsave("figures2/bray_curtis_per_species.png", width=90, height = 70, units = "mm", bg="white", dpi=300)

plot_grid(p_pop_it2, plot_grid(ggplot(), bcfig, ncol=2, labels=c("", "c")), 
          labels = "auto", ncol=1, rel_heights = c(1,1.7)) 
ggsave("figures2/figure2_Routput.png", width=200, height = 150, units = "mm", bg="white", dpi=300)










# FIGURE SI cor T1 T2 ------------------------------------

rm(list = ls())
tmp = readRDS("output/corT1T2.rds")

tmp$atlas = factor(tmp$dataset_id, levels = c(5,6,17,26),
                   labels=c("Czechia", "New York", "New Zealand", "Europe"))
zscores <- readRDS("data/z_scores.rds")
zscores <- rbindlist(zscores)
zscores[, iteration:=rep(c(1:100)), .(dataset_id, timestep, species_pair)]
zscores$atlas = factor(zscores$dataset_id, levels = c(5,6,17,26),
                       labels=c("Czechia", "New York", "New Zealand", "Europe"))
zscores$timestep[zscores$timestep==3] = 2
zscores[, state:=ifelse(z_scores<(-1.96), "segregated", ifelse(z_scores>1.96, "aggregated", "neutral")), 
        .(species_pair, timestep, atlas, iteration)]
setnames(zscores, old="timestep", new="endYear")
tmp = merge(tmp, zscores, all.x=T, by=c("atlas", "iteration", "species_pair", "endYear"))

grid.col = c("segregated"="#1398E9", "neutral"="grey", "aggregated"="#E96513")
tmp$endYear2 = factor(tmp$endYear, levels = c(1,2), labels=c("T1", "T2"))

plot_grid(labels = "auto", ncol=3, rel_widths = c(1,1.1,1), 
          ggplot(tmp, aes(x=rho, group=iteration))+
            geom_density(col="grey20", lwd=0.07)+
            stat_central_tendency(geom="line", type = "median", alpha=0.01)+
            facet_grid(atlas~endYear2, scales="free_y")+
            labs(x="Co-occurrence (\U03C1)", y="Density")+
            theme(strip.text.y = element_text(angle=-90)),
          ggplot(tmp, aes(x=rho, fill=state))+  
            geom_histogram()+
            scale_fill_manual(values=grid.col)+
            labs(x="Co-occurrence (\U03C1)", y="Species pairs")+
            stat_central_tendency(geom="line", type = "median", alpha=0.01)+
            facet_grid(atlas~endYear2, scales="free_y")+
            theme(legend.position = "bottom", legend.direction = "horizontal",
                  strip.text.y = element_text(angle=-90)),
          ggplot(tmp, aes(x=rho, y=z_scores, group=iteration))+
            geom_scattermore()+
            geom_hline(yintercept = 0, lty=2)+
            labs(x="Co-occurrence (\U03C1)", y="Z-score \U03C1")+
            facet_grid(atlas~endYear2, scales="free_y")+
            theme(strip.text.y = element_text(angle=-90))
)
ggsave("figures2/cor_T1_T2.png", width = 9, height = 3.5, dpi=150)


tmpsub = tmp[treatment=="uncertainty" & iteration==1, ]
tab = tmpsub[order(atlas, endYear2), .N, .(atlas, endYear2, state)][, perc := 100 * N / sum(N),  .(atlas, endYear2)]
tab[, perc:=round(perc,2)]
tab = gt::gt(tab)
tab

knitr::kable(tab, digits = 2)











# ---------------------------------------------------------------------------#
# Trait figures --------------------------------------------------------------
# ---------------------------------------------------------------------------#


rm(list = ls())
library(mgcv)
library(clootl)
library(ape)
library(cluster)
library(phytools)
library(performance)
library(lmerMultiMember)

# Model selection & figure ---------------------------------------------------------

files = dir("output", pattern = "^model", full.names = T)
out = lapply(files, readRDS)
names(out) = sapply(out, "[[", "atlas")

mods = sapply(out, "[", "models", USE.NAMES = T)
perf <-list()
for(i in 1:length(mods)){
  tmp <- as.data.table(unlist(lapply(mods[[i]], AIC)), keep.rownames = T)
  tmp$atlas = gsub("\\.models", "", names(mods[i]))
  perf[[i]] <- tmp
}
mod_perf <- rbindlist(perf)

### anovas ###
# Czechia
tmp <- mods$`5.models`
tmp1 =rbind(
  broom::tidy(anova(tmp$m_base, tmp$m_pl)),
  broom::tidy(anova(tmp$m_base, tmp$m_tn)),
  broom::tidy(anova(tmp$m_base, tmp$m_hab)),
  broom::tidy(anova(tmp$m_base, tmp$m_occ)))
# occ is best

# NY
tmp <- mods$`6.models`
tmp2 =rbind(
  broom::tidy(anova(tmp$m_base, tmp$m_pl)),
  broom::tidy(anova(tmp$m_base, tmp$m_tn)),
  broom::tidy(anova(tmp$m_base, tmp$m_hab)),
  broom::tidy(anova(tmp$m_base, tmp$m_occ)))
# occ is best

# NZ
tmp <- mods$`17.models`
tmp3 =rbind(
  broom::tidy(anova(tmp$m_base, tmp$m_pl)),
  broom::tidy(anova(tmp$m_base, tmp$m_tn)),
  broom::tidy(anova(tmp$m_base, tmp$m_hab)),
  broom::tidy(anova(tmp$m_base, tmp$m_occ)))
# base is best

# NY
tmp <- mods$`26.models`
tmp4 =rbind(
  broom::tidy(anova(tmp$m_base, tmp$m_pl)),
  broom::tidy(anova(tmp$m_base, tmp$m_tn)),
  broom::tidy(anova(tmp$m_base, tmp$m_hab)),
  broom::tidy(anova(tmp$m_base, tmp$m_occ)))
# occ is better than base, trophic niche is on the edge

knitr::kable(rbind(unique(tmp1), unique(tmp2), unique(tmp3), unique(tmp4)), digits = 3)  

# some models are better with occ, hence we use the same for all as we do not use it for predictions etc
bm <- list("Czechia"=mods$`5.models`$m_occ, 
           "New York"=mods$`6.models`$m_occ, 
           "New Zealand"=mods$`17.models`$m_occ, 
           "Europe"=mods$`26.models`$m_occ) 


# fixed effects table
fe = lapply(bm, broom.mixed::tidy, conf.int=TRUE, effects = c("fixed", "ran_pars")) #, "ran_vals"
fe = rbindlist(fe, use.names = T, idcol = "atlas")
fe$term2 = factor(fe$term,
                  levels = c("pd_sc", "gower_sc", "habitat_same", "min_occ_sc", 
                             "occ_ratio", "occ_change_mean_sc"),
                  labels = c("phylogenetic distance", "functional dissimilarity", 
                             "habitat agreement", "min occupancy", 
                             "occupancy ratio", "\U0394 occupancy mean"))

(p1 = ggplot(fe[term %in% c("pd_sc", "gower_sc", "habitat_same", "min_occ_sc", 
                            "occ_ratio", "occ_change_mean_sc"),], 
             aes(x=estimate, y=term2, col=!(conf.high<0|conf.low>0)))+
    geom_rect(aes(ymin = c(1.5), ymax =c(2.5), xmin = -Inf, xmax = Inf), 
              fill="grey96", col=NA)+
    geom_rect(aes(ymin = c(3.5), ymax =c(4.5), xmin = -Inf, xmax = Inf), 
              fill="grey96", col=NA)+
    geom_rect(aes(ymin = c(5.5), ymax =c(6.5), xmin = -Inf, xmax = Inf), 
              fill="grey96", col=NA)+
    geom_point()+
    geom_linerange(aes(xmin = conf.low, xmax=conf.high), lwd=0.8)+
    geom_vline(xintercept = 0,lty=2)+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          legend.background = element_blank(),
          legend.box.margin = unit(c(0,0,0,0), "mm"),
          legend.key.height = unit(5, "mm"),
          axis.ticks.y=element_blank(),
          legend.margin = margin(0,0,0,0,"mm"),
          legend.text = element_text(margin=margin(0,0,0,0,"mm")), 
          legend.box.spacing = margin(0,0,0,0,"mm")
          #legend.position.inside = c(0.3,0.12)
    )+
    labs(y="", x="Estimate")+
    scale_color_scico_d("95% CIs overlap 0", 
                        palette=scico_palette_names(categorical = FALSE)[25], begin=0.3, end=0.7, direction = -1)+
    facet_grid("atlas")
)


tmp = lapply(bm, ggeffects::ggpredict)
tmp = lapply(tmp, rbindlist, use.names = T, idcol = "variable")  
tmp = rbindlist(tmp, idcol = "atlas")  

# plot only significant ones here:
fe[, sig:=ifelse(conf.high<0, "sig", ifelse(conf.low>0, "sig", "no")),]
fesig <- fe[sig!="no"]
setnames(fesig, "term", "variable")
tmp <- merge(tmp, fesig[, .(atlas, variable, sig)], all.x=T)
tmp$variable2 = factor(tmp$variable,
                       levels = c("pd_sc", "gower_sc", "habitat_same", "min_occ_sc", 
                                  "occ_ratio", "occ_change_mean_sc"),
                       labels = c("phylogenetic distance", "functional dissimilarity", 
                                  "habitat agreement", "min occupancy", 
                                  "occupancy ratio", "\U0394occupancy mean"))
p_tmp <- tmp[sig!="no"]
p_tmp <- droplevels(p_tmp)
(p2 = ggplot(p_tmp, aes(x=x, y=predicted, col=atlas, fill=atlas))+
    geom_abline(slope = c(0), lty=2)+
    geom_line(lwd=1)+
    #scale_x_continuous(breaks = c(-4: 7))+
    labs(x="model variable", y="Predicted \U0394\U03C1")+
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.1, lty=2, show.legend = F, fill="grey80")+
    facet_wrap("variable2", scales="free_x")+
    scale_color_startrek()+
    scale_fill_startrek()+
    theme(strip.text = element_text())
)

cowplot::plot_grid(p1, p2, rel_widths = c(1,1.4), labels = "auto")
ggsave("figures2/trait_fig.png", width=8.5, height=4, dpi=300)
ggsave("figures2/trait_fig.svg", width=8.5, height=4)






## tables
# r2
tb = rbindlist(lapply(bm, r2), idcol = "atlas")
tb$R2_marginal = round(tb$R2_marginal, 2)
tb$R2_conditional = round(tb$R2_conditional, 2)
tb = gt::gt(tb)
tb


fe$estimate = round(fe$estimate, 3)
fe$std.error = round(fe$std.error, 3)
fe$conf.low = round(fe$conf.low, 3)
fe$conf.high = round(fe$conf.high, 3)

knitr::kable(fe[, .(atlas, effect, term, estimate, std.error, conf.low, conf.high)], digits = 3)

tb = gt::gt(fe[, .(atlas, effect, term, estimate, std.error, conf.low, conf.high)])
tb




