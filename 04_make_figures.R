# Script: 04_make_figures.R
# Author: Melanie Tietje
# Email: tietje@fzp.czu.cz
# GitHub: @Eryops1
# Purpose: Make figures, get stats.
# Notes: groundhog will ensure the exact same R packages will be used and create
#        a library on first run, which might take a while
#         

# Load libraries ----------------------------------------------------------

rm(list=(ls())) # clear workspace
library(groundhog)
pkgs = c("data.table",
          "ggplot2",
          "cowplot",
          "ggsci",
          "ggtext",
          "ggpubr",
          "terra",
          "sf",
          "circlize", 
          "compositions",
          "psych")
groundhog.library(pkgs, "2025-10-25")

theme_set(theme_classic(base_size = 9)+
            theme(strip.background = element_blank())+
            theme(axis.title = ggtext::element_markdown())+
            theme(axis.title.x = ggtext::element_markdown())+
            theme(axis.title.y = ggtext::element_markdown())+
            theme(legend.title = ggtext::element_markdown()))

figure_folder <- "figures_nature/"


# functions --------------------------------------------------------------

# sorts species pair names string into alphabetic order
sort_pairs_vectorized = function(x) {
  sapply(x, function(item) {
    spl = unlist(strsplit(item, split = "|", fixed = TRUE))
    paste0(sort(spl), collapse = "|")
  })
}


# function for permutation test
permdiff = function(x, group, rep = 1000, metric) {
  # x      = data.table
  # group  = column name (string) for grouping variable
  # metric = column name (string) for the metric variable
  
  stopifnot(is.data.table(x))
  stopifnot(group %in% names(x))
  stopifnot("atlas" %in% names(x))
  stopifnot(metric %in% names(x))
  
  # subset to groups with at least 3 members
  tmp = table(x[[group]])
  valid_groups = names(tmp[tmp >= 3])
  subx = x[get(group) %in% valid_groups]
  
  # permutation
  group_vec = c()
  atlas_vec = c()
  pvals = c()
  group_median = c()
  group_size = c()
  
  for (a in unique(subx$atlas)) {
    sub_atlas = subx[atlas == a]
    for (g in unique(sub_atlas[[group]])) {
      subset_group = sub_atlas[get(group) == g]
      if (nrow(subset_group) >= 3) {
        # observed statistic
        obs_median = median(subset_group[[metric]], na.rm = TRUE)
        # permutation distribution
        perm_stats = replicate(rep, {
          median(sample(sub_atlas[[metric]], nrow(subset_group)), na.rm = TRUE)
        })
        # two-sided p-value
        perm_pvalue = mean(abs(perm_stats - median(sub_atlas[[metric]], na.rm=TRUE)) >= 
                             abs(obs_median - median(sub_atlas[[metric]], na.rm=TRUE)))
        
        # collect
        group_vec = c(group_vec, g)
        atlas_vec = c(atlas_vec, a)
        pvals = c(pvals, perm_pvalue)
        group_median = c(group_median, median(sub_atlas[[metric]], na.rm = TRUE))
        group_size = c(group_size, nrow(subset_group))
      }
    }
  }
  
  # assemble results
  stat_res = data.table(
    group = group_vec,
    group_median = group_median,
    atlas = atlas_vec,
    p_value = pvals,
    n = group_size
  )
  
  # safer version for data.table
  stat_res[, significance := fifelse(.SD$p_value < 0.001, "***",
                                     fifelse(.SD$p_value < 0.01, "**",
                                             fifelse(.SD$p_value < 0.05, "*", "")))]
  
  
  # add y-position for plotting
  y_pos = numeric(nrow(stat_res))
  for (i in seq_len(nrow(stat_res))) {
    ord = stat_res$group[i]
    atl = stat_res$atlas[i]
    sub_data = subx[get(group) == ord & atlas == atl]
    if (nrow(sub_data) > 0) {
      y_pos[i] = max(sub_data[[metric]], na.rm = TRUE) + 0.02
    } else {
      y_pos[i] = NA_real_
    }
  }
  stat_res$y = y_pos
  return(list(data = subx, stats = stat_res))
}
set.seed(123)





# Read data ----------------------------------------------------------------

# co-occurrence change data per species pair
fin = readRDS("data/processed_spass.rds")

## add atlas names
atlas_names = c("26" = "Europe",
                 "5"  = "Czechia",
                 "6"  = "New York",
                 "17" = "New Zealand")
fin$atlas <- atlas_names[as.character(fin$dataset_id)]


# atlas level maps data (total richness per site, shapes)
dataset_id = unique(fin$dataset_id)
grid = lapply(paste0("data/all_scales_atlas_", dataset_id, ".gpkg"), vect)
grid = vect(grid)

dat = readRDS("data/sum_mean_psi_for_maps.rds")

shape_names = paste0("data/", c("europe", "czechia", "nz", "new_york_state"), ".gpkg")
shapes = lapply(shape_names, function(x){st_as_sf(vect(x))})
names(shapes) = c("eu", "cz", "nz", "ny")

# species level changes in range size (occupancy probability and presence grids)
range_change = readRDS("data/range_change.rds")

# atlas level mantel test results
mantel = readRDS("data/mantel_results.rds")

# species level bray curtis change
bc = readRDS("data/braycurtis.rds")






# Fig 1 maps ------------------------------------------------------------------

## merge occupancy data into the spatial grids
grid_sub_eu = grid[grid$datasetID==26,]
grid_sub_eu = grid_sub_eu[grid_sub_eu$scalingID==1,]
grid_sub_eu = grid_sub_eu[grid_sub_eu$siteID %in% dat$siteID[dat$datasetID==26]]
m = unique(dat[datasetID==26, .(siteID, sum.psi_siteID)])
grid_sub_eu = merge(grid_sub_eu, m, all.x=TRUE)

grid_sub_cz = grid[grid$datasetID==5,]
grid_sub_cz = grid_sub_cz[grid_sub_cz$scalingID==1,]
grid_sub_cz = grid_sub_cz[grid_sub_cz$siteID %in% dat$siteID[dat$datasetID==5]]
m = unique(dat[datasetID==5, .(siteID, sum.psi_siteID)])
grid_sub_cz = merge(grid_sub_cz, m, all.x=TRUE)

grid_sub_nz = grid[grid$datasetID==17,]
grid_sub_nz = grid_sub_nz[grid_sub_nz$scalingID==1,]
grid_sub_nz = grid_sub_nz[grid_sub_nz$siteID %in% dat$siteID[dat$datasetID==17]]
m = unique(dat[datasetID==17, .(siteID, sum.psi_siteID)])
grid_sub_nz = merge(grid_sub_nz, m, all.x=TRUE)

grid_sub_ny = grid[grid$datasetID==6,]
grid_sub_ny = grid_sub_ny[grid_sub_ny$scalingID==2,]
grid_sub_ny = grid_sub_ny[grid_sub_ny$siteID %in% dat$siteID[dat$datasetID==6]]
m = unique(dat[datasetID==6, .(siteID, sum.psi_siteID)])
grid_sub_ny = merge(grid_sub_ny, m, all.x=TRUE)

# map theme set
map_theme = theme_classic(base_size = 8)+
            theme(strip.background = element_blank())+
            theme(axis.title = ggtext::element_markdown())+
            theme(legend.title = ggtext::element_markdown())+
            theme_void()+
            theme(legend.position = "inside", legend.title = element_text(size=8),
                  legend.direction = "horizontal", legend.text.align=0,
                  legend.key.height = unit(4,"mm"), legend.key.width = unit(5,"mm"),
                  legend.background = element_blank(),legend.key = element_blank(),
                  panel.background = element_blank(), legend.position.inside = c(0.2,0.9),
                  legend.title.position = "top")
lw = 0.2

# maps
(eu_map = ggplot(shapes$eu)+
    geom_sf(data=st_as_sf(grid_sub_eu), aes(fill=sum.psi_siteID), col=NA)+
    geom_sf(col="black", fill="NA", lwd=lw)+
    scale_fill_gradient("Species richness", low = 'white', high = '#5f89d8',
                        breaks = c(min(grid_sub_eu$sum.psi_siteID),max(grid_sub_eu$sum.psi_siteID)),
                        labels = c(round(min(grid_sub_eu$sum.psi_siteID)),round(max(grid_sub_eu$sum.psi_siteID))))+
  map_theme)

(cz_map = ggplot(shapes$cz)+
    geom_sf(data=st_as_sf(grid_sub_cz), aes(fill=sum.psi_siteID), col=NA)+
    geom_sf(col="black", fill="NA", lwd=lw)+
    scale_fill_gradient("Species richness", low = 'white', high = '#c90f04',
                        breaks = c(min(grid_sub_cz$sum.psi_siteID),max(grid_sub_cz$sum.psi_siteID)),
                        labels = c(round(min(grid_sub_cz$sum.psi_siteID)),round(max(grid_sub_cz$sum.psi_siteID))))+
    map_theme)

(nz_map = ggplot(shapes$nz)+
    geom_sf(data=st_as_sf(grid_sub_nz), aes(fill=log(sum.psi_siteID)), col=NA)+
    geom_sf(col="black", fill="NA", lwd=lw)+
    scale_fill_gradient("Species richness", low = 'white', high = '#fecd04',##fecd04',#ebbe04
                        breaks = c(min(log(grid_sub_nz$sum.psi_siteID)),max(log(grid_sub_nz$sum.psi_siteID))),
                        labels = c(round(min(grid_sub_nz$sum.psi_siteID)),round(max(grid_sub_nz$sum.psi_siteID))))+
    map_theme)

(ny_map = ggplot(shapes$ny)+
    geom_sf(data=st_as_sf(grid_sub_ny), aes(fill=sum.psi_siteID), col=NA)+
    geom_sf(col="black", fill="NA", lwd=lw)+
    scale_fill_gradient("Species richness", low = 'white', high = '#85bb07',
                        breaks = c(min(grid_sub_ny$sum.psi_siteID),max(grid_sub_ny$sum.psi_siteID)),
                        labels = c(round(min(grid_sub_ny$sum.psi_siteID)),round(max(grid_sub_ny$sum.psi_siteID))))+
    map_theme)

# assemble
plot_grid(plot_grid(eu_map+theme(plot.background = element_rect(), title = element_text(size=8))+ggtitle("Europe"), 
                    nz_map+theme(plot.background = element_rect(), title = element_text(size=8))+ggtitle("New Zealand"), 
                    rel_widths = c(3,2.63)),  
          plot_grid(cz_map+theme(plot.background = element_rect(), title = element_text(size=8))+ggtitle("Czechia"), 
                    ny_map+theme(plot.background = element_rect(), title = element_text(size=8))+ggtitle("New York State"), 
                    ncol = 1, rel_heights = c(1,0.9)), 
          ncol=2, rel_widths = c(3,1.2))
# doing the fine tuning and arranging in inkscape
ggsave(paste0(figure_folder, "data_maps.svg"), width=7.6, height=4, bg = "white")
rm(shapes, grid_sub_cz, grid_sub_eu, grid_sub_ny, grid_sub_nz, m, map_theme)
rm(list=ls(pattern="_map"))





# Fig 1 timeline ---------------------------------------------------------------
# number of species and pairs in which years
## # test if species pair overlaps in all sampling periods & subset
fin[, over:=ifelse(all(overlap==FALSE), FALSE, TRUE), by=species_pair]
temp = unique(fin[,.(atlas, time_bin, over, species_pair)])
temp = temp[over==TRUE,] # subset to pairs that overlap

temp$years = temp$time_bin
species = temp[,unique(unlist(strsplit(species_pair, split = "|", fixed = TRUE))), by=atlas]
species = species[,length(V1), by=atlas]
tmp = data.table(table(temp$atlas[temp$time_bin=="T1"]))
setnames(tmp, old=c('V1', 'N'), new=c('atlas', 'N pairs'))

# add years for sampling periods
cz = c(c(1985:1989), c(2001:2003), c(2014:2017))
cz = data.table(year = cz, atlas = rep("Czechia", length(cz))) 
ny = c(c(1980:1985), c(2000:2005)) # 1980-1985, 2000-2005
ny = data.table(year = ny, atlas = rep("New York", length(ny))) 
nz = c(c(1969:1979), c(1999:2004)) # 1969-1979, 1999-2004
nz = data.table(year = nz, atlas = rep("New Zealand", length(nz))) 
eu = c(c(1972:1995), c(2013:2017)) # 1972-1995, 2013-2017
eu = data.table(year = eu, atlas = rep("Europe", length(eu))) 
temp = rbind(cz,eu,ny,nz)
temp = merge(temp, tmp)
temp = merge(temp, species, by='atlas', all.x=TRUE)

ggplot(temp, aes(x=year, y=atlas, col=atlas, size=2))+
  geom_point(pch=15)+
  scale_color_startrek(guide='none')+
  scale_size_continuous(guide='none')+
  labs(y="")+
  theme(panel.background = element_blank())
# separate dots here but saving it small enough will connect them just fine
ggsave(paste0(figure_folder, "data_timeline.svg"), height=1.5, width=2.4, bg = NULL)

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


rm(temp, cz, ny, nz, eu)




# stats and circlize plot -----------------------------------------------------

stat = fin[over==TRUE,] # subset to pairs that overlap

# add status based on Z-score (more or less aggregated than under null distribution)
stat[,status_T1:=cut(cor_ses[time_bin=="T1"], 
                     breaks=c(-Inf,-1.96,1.96,Inf), 
                     labels=c("segregated", "neutral", "aggregated")), 
     by=.(atlas, species_pair, scaleID)]
stat[,status_Tlast:=cut(cor_ses[time_bin==tail(time_bin,1)], 
                        breaks=c(-Inf,-1.96,1.96,Inf), 
                        labels=c("segregated", "neutral", "aggregated")), 
     by=.(atlas, species_pair, scaleID)]
stat = unique(stat[,.(species_pair, atlas, status_T1, status_Tlast)])

# proportions 
d1 = as.data.table(table(stat$status_T1, stat$atlas))
d1$sampling_period = "first"
d2 = as.data.table(table(stat$status_Tlast, stat$atlas))
d2$sampling_period = "last"
dt = rbind(d1,d2)

# add percent
dt[, prop:=N/sum(N), by=.(sampling_period, V2)]
dt[, prop_per_group:=mean(prop), by=.(sampling_period, V1)]
dt[, prop_per_group_range:=max(prop)-min(prop), by=.(sampling_period, V1)]
dt[, prop_change:=prop[sampling_period=="last"]-prop[sampling_period=="first"], by=.(V1, V2)]
setnames(dt, c('V1', 'V2'), c('co-occurrence_status', 'atlas'))
# CZ: 1985-1989, 2001-2003, 2014-2017
# NY: 1980-1985, 2000-2005
# EU: 1972-1995, 2013-2017
# NZ: 1969-1979, 1999-2004
sp = data.table(atlas = rep(c('Czechia', 'New York', 'Europe', 'New Zealand'), each=2), 
                sampling_period = rep(c('first', 'last'), 4), 
                years = c('1985-1989', '2014-2017', '1980-1985', '2000-2005', 
                          '1972-1995', '2013-2017', '1969-1979', '1999-2004'))
dt = merge(dt,sp)

fwrite(dt[, .(atlas, `co-occurrence_status`, sampling_period, years, N, 
              round(prop,3), round(prop_per_group,3), round(prop_per_group_range,3), round(prop_change,3))], 
       "output/status_stats1.csv")
rm(d1,d2,sp)





## circlize plot ---- 

tmp = unique(stat[,.(species_pair, atlas, status_T1, status_Tlast)])
tmp[,transition:=paste(status_T1, status_Tlast, sep = "-->")]

tmp2 = as.data.table(table(tmp$transition, tmp$atlas))
tmp2[, status_T1:=gsub("-->.*", "", V1)]
tmp2[, status_T2:=gsub(".*-->", "", V1)]
tmp2[, change_status:=status_T1!=status_T2]
tmp2[, sum_per_start:=sum(N), by=.(V2, status_T1)]
tmp2[, prop_per_trans:=N/sum_per_start, by=.(V2, status_T1)]

tmp3 = tmp2[, sum(N), by=.(change_status)]
tmp3[, prop:=round(V1/sum(V1), 3)]
fwrite(tmp3, "output/transition_stats1.csv")
tmp3 = tmp2[, sum(N), by=.(change_status, V2)]
tmp3[, prop:=round(V1/sum(V1), 3), by=V2]
fwrite(tmp3, "output/transition_stats2.csv")
fwrite(tmp2, "output/transition_stats3.csv")

# atli
atli = sort(unique(tmp$atlas))
grid.col = c("segregated"="#1398E9", "neutral"="grey", "aggregated"="#E96513")

png(paste0(figure_folder, "circlize.png"), width =11, height=11.4, units = "cm", bg = "transparent", res = 300)
par(mfrow=c(2,2))
for(a in 1:length(atli)){
  tmp_cz = as.data.table(table(tmp$transition[tmp$atlas==atli[a]]))
  tmp_cz[, from:=gsub("-->.*", "", V1)]
  tmp_cz[, to:=gsub(".*-->", "", V1)]
  setnames(tmp_cz, 'N', 'value')
  tmp_cz = tmp_cz[,.(from, to, value)]
  chordDiagram(tmp_cz,
               directional = 1, direction.type = c("arrows"),
               link.arr.type= "big.arrow",
               #link.arr.length = ifelse(link.arr.type == "big.arrow", 0.02, 0.4),
               link.arr.length = 0.1,
               self.link = 1,
               link.largest.ontop=F, grid.col = grid.col, #link.border = "grey40", link.lwd = 1,
               transparency=0.2,
               link.target.prop = TRUE,
               #target.prop.height=0, 
               annotationTrack = c("grid"),
               annotationTrackHeight = c(grid = 0.06, name=0.13),
               order = c("aggregated", "neutral", "segregated"))
  # transition percentages per category
  cats = unique(tmp_cz$from)
  cat_props = round(tmp_cz$value/sum(tmp_cz$value)*100, 1)
  names(cat_props) = tmp_cz$from
  
  # sum percentages per category (total across transitions)
  cat_totals = sapply(cats, function(cat) {
    sum(cat_props[grepl(cat, names(cat_props))], na.rm = TRUE)
  })
  circos.track(track.index = 1, panel.fun = function(x, y) {
    sector = get.cell.meta.data("sector.index")
    
    x_pos = mean(get.cell.meta.data("xlim"))
    y_pos = get.cell.meta.data("ylim")[2] + 1.2  # place above the sector
    if (sector %in% names(cat_totals)) {
      label = paste0(sector, " (", format(cat_totals[sector], nsmall = 1), "%)")
      circos.text(x = x_pos, y = y_pos,
                  labels = label,
                  facing = "bending.inside", niceFacing = TRUE, cex = 0.9, col = "black")
    }
  },
  , bg.border = NA)
  title(main = atli[a], cex.main=1.1, font.main=1)
}
dev.off()









# correlation histograms per sampling period -------------------------------

# remove middle sampling period for CZ & rename last
fig_hist = fin[!(time_bin=="T2"& atlas=="Czechia"),] 
fig_hist$time_bin[fig_hist$time_bin=="T3"] = "T2"
fig_hist$x_bin = cut(
  fig_hist$cor_ses,
  breaks = c(-Inf, -1.96, 1.96, Inf),
  labels = c("Seg (< -1.96)", "Neutral (-1.96 to 1.96)", "Agg (> 1.96)"),
  ordered_result = TRUE   # ensures ordering
)
cor_per_time <- ggplot(fig_hist[over==TRUE,], aes(x=cor_obs))+
  geom_histogram(binwidth = 0.1, aes(fill=x_bin), bins=160)+ # aes(y = ..density..), , fill="grey80"
  facet_grid(time_bin~atlas, scales="free_y")+
  geom_vline(xintercept = 0, lty=1, col="grey40", lwd=0.3)+
  stat_central_tendency(type="median", col="black", lty=2, lwd=0.3)+
  scale_fill_manual(name = "", values = c("#1398E9", "grey", "#E96513"), labels=c("segregated", "neutral", "aggregated"))+
  scale_color_startrek(guide="none")+
  labs(y="Count", x="Co-occurrence (Spearman's \U03C1)")+
  theme(legend.position = "bottom", strip.text.y = element_text(angle=0), strip.background.y = element_rect(fill = "grey95", linewidth = 0))

ggsave(cor_per_time, paste0(figure_folder, "cor_obs_histogram.png"), width=6, height = 2.5, bg="white")

# ggplot(fig_hist[over==TRUE,], aes(x=cor_obs))+
#   geom_histogram(binwidth = 0.1, aes(fill=atlas), bins=160)+ # aes(y = ..density..), , fill="grey80"
#   facet_grid(time_bin~atlas, scales="free_y")+
#   geom_vline(xintercept = 0, lty=1, col="grey40", lwd=0.3)+
#   stat_central_tendency(type="median", col="black", lty=2, lwd=0.3)+
#   scale_fill_startrek(guide="none")+
#   scale_color_startrek(guide="none")+
#   labs(y="Count", x="Co-occurrence (Spearman's \U03C1)")+
#   theme(legend.position = "bottom", strip.text.y = element_text(angle=0), strip.background.y = element_rect(fill = "grey95", linewidth = 0))
# ggsave(paste0(figure_folder, "cor_obs_histogram_TALK.png"), width=6, height = 2.5, bg="white")

### stats ----
psych::describeBy(data=fig_hist, cor_obs ~ atlas + time_bin, mat=TRUE, digits=2)
kruskal.test(fig_hist$cor_obs[fig_hist$dataset_id==5], fig_hist$time_bin[fig_hist$dataset_id==5])
kruskal.test(fig_hist$cor_obs[fig_hist$dataset_id==6], fig_hist$time_bin[fig_hist$dataset_id==6])
kruskal.test(fig_hist$cor_obs[fig_hist$dataset_id==17], fig_hist$time_bin[fig_hist$dataset_id==17])
kruskal.test(fig_hist$cor_obs[fig_hist$dataset_id==26], fig_hist$time_bin[fig_hist$dataset_id==26])






# Z-scores histogram per sampling period ------------------------------------

fig_hist$x_bin = cut(
  fig_hist$cor_ses,
  breaks = c(-Inf, -1.96, 1.96, Inf),
  labels = c("Seg (< -1.96)", "Neutral (-1.96 to 1.96)", "Agg (> 1.96)"),
  ordered_result = TRUE   # ensures ordering
)
zhist = ggplot(fig_hist[over==TRUE,], aes(x=cor_ses))+
  geom_histogram(aes(fill=x_bin), bins=40)+ # aes(y = ..density..), , fill="grey80"
  facet_grid(time_bin~atlas, scales="free_y")+
  geom_vline(xintercept = 0, lty=1, col="grey40", lwd=0.3)+
  stat_central_tendency(type="median", col="black", lty=2, lwd=0.3)+
  #scale_fill_startrek(guide="none")+
  #grid.col = c("segregated"="#1398E9", "neutral"="grey", "aggregated"="#E96513")
  scale_fill_manual(name = "", values = c("#1398E9", "grey", "#E96513"), labels=c("segregated", "neutral", "aggregated"))+
  scale_color_startrek(guide="none")+
  labs(y="Count", x="Co-occurrence (_Z_-score Spearman's \U03C1)")+
  theme(legend.position = "bottom", 
        strip.text.y = element_text(angle=0), 
        strip.background.y = element_rect(fill = "grey95", linewidth = 0),
        panel.spacing.x = unit(0.3, "cm"))

### stats -----
psych::describeBy(data=fin, cor_ses ~ atlas + time_bin, mat=TRUE, digits=2)








## raw vs Z-scores ----
corscat = ggplot(fig_hist, aes(x=cor_obs, y=cor_ses, col=atlas))+
  geom_point(alpha=0.01)+
  stat_cor(method = "spearman", label.x = 0.2, label.y = -25, size=2.5, 
           label.sep = "\n", p.accuracy = 0.001, cor.coef.name = "rho", col="black")+
  facet_grid(time_bin~atlas)+
  scale_color_startrek(guide="none")+
  labs(y="*Z*-score Spearman's \U03C1", x="Spearman's \U03C1")+
  theme(legend.position = "bottom", 
        strip.text.y = element_text(angle=0), 
        strip.background.y = element_rect(fill = "grey95", linewidth = 0),
        panel.spacing.x = unit(0.3, "cm"))

plot_grid(cor_per_time, zhist+theme(legend.position = "none"), corscat, 
          ncol=1, labels = "auto", label_fontface = "bold", 
          label_size = 10, rel_heights = c(1,0.85,1))
ggsave(paste0(figure_folder, "figS1.png"), width=6, height = 6.5, bg="white")

rm(zhist, corscat)



































# PAIR CHANGE  ------------------------------------------------------------------

fin[, cor_obs_T1:=cor_obs[time_bin=="T1"], by=.(atlas, scaleID, species_pair)]
chan = unique(fin[,.(over, species_pair, atlas, delta_cor_obs, 
                     delta_c_obs, delta_c_obs_z, delta_cor_obs_z, cor_obs_T1)])
chan = chan[over==TRUE,]



## get change stats for plot
chan[,delta_cor_obs_median:=median(delta_cor_obs, na.rm=TRUE), by=atlas]
chan[,delta_cor_obs_mean:=mean(delta_cor_obs, na.rm=TRUE), by=atlas]
statdat = unique(chan[,.(delta_cor_obs_median, delta_cor_obs_mean,atlas)])
statdat$delta_cor_obs_median = round(statdat$delta_cor_obs_median,3)
statdat$delta_cor_obs_mean = round(statdat$delta_cor_obs_mean,3)
statdat$label_med = paste0('median:\n', statdat$delta_cor_obs_median)
statdat$label_m = paste0('mean:\n', statdat$delta_cor_obs_mean)


## plot
(fig3a = ggplot(chan, aes(x=delta_cor_obs, col=atlas, fill=atlas))+
  geom_density()+
  facet_wrap(~atlas, ncol=4)+
  labs(y="Density", x="Co-occurrence change\n(\U0394 Spearman's \U03C1)")+
  stat_central_tendency(type="median", aes(col=atlas), lty=2)+
  scale_color_startrek(guide='none')+
  scale_fill_startrek(guide='none', alpha=0.5)+
  geom_text(data=statdat, aes(label=label_med), y=2, x=0.9, col="black", size=2.3, hjust=1))
ggsave(paste0(figure_folder, "delta_cor_obs_density.png"), width=6, height = 1.7, bg="white", dpi=300)

### stats -----
# by atlas
psych::describeBy(data=chan, delta_cor_obs ~ atlas, mat=TRUE, digits=2, IQR=TRUE, quant=c(.25,.75))
tapply(chan$delta_cor_obs, chan$atlas, function(x){shapiro.test(sample(x, 5000, replace = FALSE))})
# all data
psych::describe(chan$delta_cor_obs, IQR=TRUE, quant=c(.25,.75))


# # V2, for poster
# ggplot(chan, aes(x=delta_cor_obs, col=atlas, fill=atlas))+
#   geom_density()+
#   facet_wrap(~atlas, ncol=2)+
#   #scale_y_sqrt()+
#   labs(y="density", x="Co-occurrence change\n(\U0394 Spearman's \U03C1)")+
#   stat_central_tendency(type="mean", aes(col=atlas), lty=2)+
#   scale_color_startrek(guide='none')+
#   scale_fill_startrek(guide='none', alpha=0.5)+
#   geom_text(data=statdat, aes(label=label_m), y=2, x=0.5, col="black", size=2.5)
# ggsave(paste0(figure_folder, "delta_cor_obs_density_2x2.png"), width=3, height = 3, bg="white", dpi=600)








## Z-score change, SI ----

## get change stats for plot
chan[,delta_cor_obs_z_median:=median(delta_cor_obs_z, na.rm=TRUE), by=atlas]
chan[,delta_cor_obs_z_mean:=mean(delta_cor_obs_z, na.rm=TRUE), by=atlas]
statdat = unique(chan[,.(delta_cor_obs_z_median, delta_cor_obs_z_mean,atlas)])
statdat$delta_cor_obs_z_median = round(statdat$delta_cor_obs_z_median,3)
statdat$delta_cor_obs_z_mean = round(statdat$delta_cor_obs_z_mean,3)
statdat$label_med = paste0('median:\n', statdat$delta_cor_obs_z_median)
statdat$label_m = paste0('mean:\n', statdat$delta_cor_obs_z_mean)

ggplot(chan, aes(x=delta_cor_obs_z, col=atlas, fill=atlas))+
    geom_density()+
    facet_wrap(~atlas, ncol=4)+
    labs(y="Density", x="Co-occurrence change (*Z*-score \U0394 Spearman's \U03C1)")+
    stat_central_tendency(type="mean", aes(col=atlas), lty=2)+
    scale_color_startrek(guide='none')+
    scale_fill_startrek(guide='none', alpha=0.5)+
    geom_text(data=statdat, aes(label=label_med), y=0.04, x=30, col="black", size=2.5, hjust=1)
ggsave(paste0(figure_folder,  "delta_cor_obs_z_density.png"), width=6, height = 1.7, bg="white", dpi=300)


### stats ----
psych::describeBy(data=chan, delta_cor_obs_z ~ atlas, mat=TRUE, digits=2)














  
  
  








# MANTEL test results --------------------------------------------------------

mantel_PD_res = mantel$mantel_PD_res
mantel_FD_res = mantel$mantel_FD_res


### PD ----
perms = as.data.table(sapply(mantel_PD_res, "[[", "perm"))
perms = melt(perms, measure.vars = names(perms), value.name = "PD_perm", variable.name = "datasetID")
stat = data.table(statistic = sapply(mantel_PD_res, "[[", "statistic"),
                   p_value = sapply(mantel_PD_res, "[[", "signif"),
                      datasetID = names(sapply(mantel_PD_res, "[[", "statistic")))
perms = merge(perms, stat)
perms$atlas = perms$datasetID
perms$atlas <- atlas_names[as.character(perms$atlas)]
perms$metric = "PD"

pd_tab = unique(perms[, .(datasetID, statistic, p_value, atlas)])

keep = perms

### FD -----
perms = as.data.table(sapply(mantel_FD_res, "[[", "perm"))
perms = melt(perms, measure.vars = names(perms), value.name = "FD_perm", variable.name = "datasetID")
stat = data.table(statistic = sapply(mantel_FD_res, "[[", "statistic"),
                   p_value = sapply(mantel_PD_res, "[[", "signif"),
                   datasetID = names(sapply(mantel_FD_res, "[[", "statistic")))
perms = merge(perms, stat)
perms$atlas = perms$datasetID

perms$atlas <- atlas_names[as.character(perms$atlas)]
perms$metric = "FD"

# combine PD and FD
keep = rbind(keep, perms, use.names=FALSE)
fd_tab = unique(perms[, .(datasetID, statistic, p_value, atlas)])



# table
mantel_results_pd = data.table(
  atlas_ID = names(mantel_PD_res),
  Matrix1 = rep("spatial association", 4),
  Matrix2 = rep("phylogenetic distance", 4),
  Mantel_r = round(sapply(mantel_PD_res, function(x) x$statistic),3),
  p_value = sapply(mantel_PD_res, function(x) x$signif),
  Permutations = sapply(mantel_PD_res, function(x) x$permutations)
)
mantel_results_fd = data.table(
  atlas_ID = names(mantel_FD_res),
  Matrix1 = rep("spatial association", 4),
  Matrix2 = rep("functional distance", 4),
  Mantel_r = round(sapply(mantel_FD_res, function(x) x$statistic), 3),
  p_value = sapply(mantel_FD_res, function(x) x$signif),
  Permutations = sapply(mantel_FD_res, function(x) x$permutations)
)
tab = rbind(mantel_results_pd, mantel_results_fd)
tab$atlas = as.character(tab$atlas_ID)
tab$atlas <- atlas_names[as.character(tab$atlas)]
tab$significance = ifelse(tab$p_value < 0.001, "***",
                          ifelse(tab$p_value < 0.01, "**",
                                 ifelse(tab$p_value < 0.05, "*", "")))

#write.csv(tab, "output/mantel_results.csv")

# add positions for asterisks placement
tab$y_pos = rep(c(100,140,140,80),2)
tab$x_pos = tab$Mantel_r
tab$metric = ifelse(tab$Matrix2 == "phylogenetic distance", "PD", "FD")


mantel_fig <- ggplot(keep, aes(x=PD_perm, fill=atlas))+
  geom_histogram(alpha=0.7)+
  geom_vline(aes(xintercept = statistic), lty=2, col="black")+
  scale_fill_startrek(guide="none")+
  scale_color_startrek(guide="none")+
  stat_central_tendency(type="mean", aes(col=atlas), lty=1)+
  labs(x="Mantel test score (\U03C1)", y="Count")+
  facet_grid(atlas~metric, scales="free_y")+
  geom_text(data=tab, aes(x=Mantel_r+0.01, y=y_pos, label = significance),
            inherit.aes = FALSE,
            hjust = 0,
            size = 3)+
  theme(strip.text = element_text(size=6))
# tab2 <- tab[, .(atlas, Mantel_r, p_value, significance, metric)]
# setnames(tab2, old=c("atlas", "metric", "Mantel_r", "p_value"), 
#          new=c("Atlas", "Metric", "Mantel r", "P-value"),skip_absent = T)
# mantel_fig + wrap_table(tab2[, .(Atlas, `Mantel r`, `P-value`, Metric)]) + 
#   plot_layout(ncol=2)

ggsave(paste0(figure_folder,  "mantel_histograms.png"), width=140, height = 72, 
       unit="mm", bg="white", dpi=300)

rm(list=ls(pattern="mantel"))
























# SPECIES LEVEL -------------------------------------------------------------

# get average values for individual species from all pairs

indiv = list()
indiv2 = list()
at = sort(unique(chan$atlas))
for(a in 1:length(at)){
  tmp = chan[atlas==at[a],]
  tmp$sp1 = gsub("\\|.*", "", tmp$species_pair)
  tmp$sp2 = gsub(".*\\|", "", tmp$species_pair)
  indiv_int = data.table(species = sort(unique(c(tmp$sp1, tmp$sp2))), 
                         delta_cor_obs_median = NA,
                         delta_c_obs_median = NA,
                         atlas = at[a],
                         n_pairs = NA)
  for(i in 1:nrow(indiv_int)){
    id = grep(indiv_int$species[i], tmp$species_pair)
    indiv_int$delta_cor_obs_median[i] = median(tmp$delta_cor_obs[id], na.rm=TRUE)
    indiv_int$delta_cor_obs_sum[i] = sum(tmp$delta_cor_obs[id], na.rm=TRUE)
    indiv_int$delta_cor_obs_mean[i] = mean(tmp$delta_cor_obs[id], na.rm=TRUE)
    indiv_int$delta_cor_obs_sd[i] = sd(tmp$delta_cor_obs[id], na.rm=TRUE)
    indiv_int$cor_obs_T1_median[i] = median(tmp$cor_obs_T1[id], na.rm=TRUE)
    indiv_int$delta_c_obs_median[i] = median(tmp$delta_c_obs[id], na.rm=TRUE)
    indiv_int$n_pairs[i] = length(id)
    if(i==1){
      indiv_int2 = data.table(delta_cor_obs = tmp$delta_cor_obs[id],
                              species = indiv_int$species[i],
                              atlas = at[a],
                              species_pair = tmp$species_pair[id])
    }else{tmp2 = data.table(delta_cor_obs = tmp$delta_cor_obs[id],
                           species = indiv_int$species[i],
                           atlas = at[a],
                           species_pair = tmp$species_pair[id])
      indiv_int2 = rbind(indiv_int2, tmp2)}

  }
  indiv[[a]] = indiv_int
  indiv2[[a]] = indiv_int2
  cat(a, i,  "\r")
}
indiv = rbindlist(indiv)
indiv2 = rbindlist(indiv2) 
# same as indiv, but keeps all values per species not just the averages
rm(indiv_int, indiv_int2)




## bray curtis dissimilarity results -----

bc$atlas <- atlas_names[as.character(bc$atlas)]
# add correlation values 
bc = merge(bc, indiv[,.(species, atlas, delta_cor_obs_mean)])
ggplot(bc, aes(x=bray_curtis, y=delta_cor_obs_mean))+geom_point()+geom_smooth(method="lm")+scale_x_continuous(trans = "sqrt")
cor.test(bc$bray_curtis, bc$delta_cor_obs_mean, method = "s")

(bcfig = ggplot(bc, aes(x=bray_curtis))+
    geom_histogram(fill="grey80")+
    facet_wrap('atlas')+
    scale_fill_startrek(guide="none")+
    scale_color_startrek(guide="none")+
    stat_central_tendency(type="median", aes(col=atlas), lty=2)+
    labs(x='Bray-Curtis between sampling periods, per species', y="Count"))
ggsave(paste0(figure_folder,  "bray_curtis_per_species.png"), width=90, height = 70, units = "mm", bg="white", dpi=300)


# add to indiv object
indiv = merge(indiv, bc, all.x=TRUE)










## TRAITS (AVONET) -------------------------------------------------------------
avonet = fread("data/AVONET1_BirdLife.csv")
names(avonet)
avonet$species = gsub(" ", "_", avonet$Species1)
indiv = merge(indiv, 
              avonet[,.SD, .SDcols = !c('Mass.Refs.Other', 'Reference.species', 'Traits.inferred')],
              by="species", all.x=TRUE)

### Primary.Lifestyle -----
stat_res = permdiff(indiv, group = "Primary.Lifestyle", rep=10000, metric="delta_cor_obs_median")
fwrite(stat_res$stats, "output/lifestyle_permutation_stats.csv")

ggplot(stat_res$data, aes(y = Primary.Lifestyle, x = delta_cor_obs_median, fill = atlas)) +
  geom_vline(data=stat_res$stats, aes(xintercept = group_median), col = "grey20", lty=2) +
  #geom_boxplot(varwidth = F, outlier.alpha = 0.2) +
  geom_boxplot(varwidth = T, outlier.alpha = 0.2) +
  scale_fill_startrek() +
  labs(y = "", x = "Species median \u0394 Spearman's \u03C1") +
  theme(legend.position = "none") +
  facet_grid(~atlas, scales="free") +
  geom_text(data = stat_res$stats, aes(x = y, y = group, label = significance),
            hjust = 0, size = 3)
# geom_text(data = stat_res$stats, aes(x = y, y = group, label = n),
#           hjust = 0, size = 3)
ggsave(paste0(figure_folder,  "delta_cor_obs_per_Primary.Lifestyle.png"), width=180, height = 40, unit="mm", bg="white", dpi=300)


# SI
### Taxonomy -----
stat_res = permdiff(indiv, group="Order1", rep=1000, metric="delta_cor_obs_median")
fwrite(stat_res$stats, "output/taxonomy_permutation_stats.csv")

taxplot = ggplot(stat_res$data, aes(y = Order1, x = delta_cor_obs_median, fill = atlas)) +
  geom_vline(data=stat_res$stats, aes(xintercept = group_median), col = "grey") +
  geom_boxplot(varwidth = F, outlier.alpha = 0.2) +
  scale_fill_startrek() +
  labs(y = "", x = "Species median \u0394 Spearman's \u03C1") +
  theme(legend.position = "none") +
  facet_wrap(~atlas, ncol = 4, scales="free_x") +
  geom_text(data = stat_res$stats,
            aes(x = y, y = group, label = significance),
            inherit.aes = FALSE,
            hjust = 0,
            size = 4)


### Habitat -----
stat_res = permdiff(indiv, group="Habitat", rep=1000, metric="delta_cor_obs_median")
fwrite(stat_res$stats, "output/habitat_permutation_stats.csv")

habplot = ggplot(stat_res$data, aes(y = Habitat, x = delta_cor_obs_median, fill = atlas)) +
  geom_vline(data=stat_res$stats, aes(xintercept = group_median), col = "grey") +
  geom_boxplot(varwidth = F, outlier.alpha = 0.2) +
  scale_fill_startrek() +
  labs(y = "", x = "Species median \u0394 Spearman's \u03C1") +
  theme(legend.position = "none") +
  facet_wrap(~atlas, ncol = 4, scales="free_x") +
  geom_text(data = stat_res$stats,
            aes(x = y, y = group, label = significance),
            inherit.aes = FALSE,
            hjust = 0,
            size = 4)

table(indiv$Habitat, indiv$atlas)





### Trophic Niche ----
stat_res = permdiff(indiv, group = "Trophic.Niche", rep=10000, metric="delta_cor_obs_median")
fwrite(stat_res$stats, "output/trophicniche_permutation_stats.csv")

trophplot = ggplot(stat_res$data, aes(y = Trophic.Niche, x = delta_cor_obs_median, fill = atlas)) +
  geom_vline(data=stat_res$stats, aes(xintercept = group_median), col = "grey") +
  geom_boxplot(varwidth = F, outlier.alpha = 0.2) +
  scale_fill_startrek() +
  labs(y = "", x = "Species median \u0394 Spearman's \u03C1") +
  theme(legend.position = "none") +
  facet_grid(~atlas, scales="free") +
  geom_text(data = stat_res$stats, aes(x = y, y = group, label = significance),
            hjust = 0, size = 4)


### PLOT
set_null_device("png") # this prevents the space after unicode to dissapear
plot_grid(taxplot, habplot, trophplot, ncol=1, labels = "auto",  
          label_fontface = "bold", label_size = 10, rel_heights = c(65, 40, 40))
ggsave(paste0(figure_folder,  "figS5.png"), width=200, height = 65+40+45+20, unit="mm", bg="white", dpi=300)
rm(taxplot, trophplot, habplot)













# RANGE SIZE  related figures -------------------------------------------------

## changes in range size and changes in correlation? ----

# merge range change into indiv data
range_change$atlas = range_change$datasetID
range_change$atlas <- atlas_names[as.character(range_change$atlas)]
range_change$species = range_change$scientificName
# add number of grids  / sum occupancy first time period
range_change[, pres.abs_sum_T1:=pres.abs_sum[endYear==1], by=.(atlas, species)]
range_change[, occupancy_sum_T1:=mean.psi_sum[endYear==1], by=.(atlas, species)]

indiv = merge(indiv, unique(range_change[,.(species, atlas, change_occupancy, change_pres.abs, 
                                            pres.abs_sum_T1, occupancy_sum_T1)]),
              by=c("atlas", "species"), 
              all.x=TRUE)

# use log ratio for change occupancy
(scatter_delta_occ <- ggplot(indiv, aes(x=as.numeric(clr(change_occupancy)), y=delta_cor_obs_median, col=atlas))+
    geom_vline(xintercept = 0, col='grey', lty=2) + geom_hline(yintercept = 0, col='grey', lty=2)+
    geom_point(alpha=0.7)+
    scale_color_startrek(guide="none")+
    geom_smooth(alpha=0.5, lwd=0.5)+
    stat_cor(method = "spearman", size=2.5, p.accuracy = 0.001, cor.coef.name = "rho",
             label.sep = "\n", col="black", label.y = 0.22)+
    facet_wrap(~atlas, ncol=4)+
    labs(x="\U0394 occupancy (centered log-ratio)", y="Median species \U0394 Spearman's \U03C1")+
    theme(legend.position = "bottom"))
ggsave(paste0(figure_folder,  "delta_occupancy_delta_correlation.png"), units="mm", width=120, height = 100, bg="white", dpi = 300)


cor.test(indiv$change_occupancy, indiv$delta_cor_obs_median, method="s")
tapply(indiv, indiv$atlas, function(x){cor.test(x$delta_cor_obs_median, as.numeric(clr(x$change_occupancy)), method='s')})





## general range size within Primary Lifestyle groups (sum occupancy all cells T1) ----

range_diff_stats = permdiff(indiv, group="Primary.Lifestyle", rep=1000, metric="occupancy_sum_T1")
fwrite(range_diff_stats$stats, "output/range_permutation_stats.csv")
(range_diff_t1 <- ggplot(range_diff_stats$data, aes(y = Primary.Lifestyle, x = occupancy_sum_T1, fill = atlas)) +
  geom_vline(data=range_diff_stats$stats, aes(xintercept = group_median), col = "grey") +
  geom_boxplot(varwidth = T, outlier.alpha = 0.2) +
  scale_fill_startrek() +
  labs(y = "", x = "Sum of occupancy probabilites in T1") +
  theme(legend.position = "none") +
  facet_grid(~atlas, scales="free")+
  geom_text(data = range_diff_stats$stats, aes(x = y, y = group, label = significance),
            hjust = 0, size = 4))
#ggsave(paste0(figure_folder,  "sum_occupancy_T1_lifestyle.png"), units="mm", width=180, height = 40, bg="white", dpi = 300)

# groups should be compared among each other and not betwen each other. permutation test
# tmp = tapply(indiv, indiv$atlas, function(x){pairwise.wilcox.test(x$occupancy_sum_T1, x$Primary.Lifestyle, p.adjust.method = "fdr")})
# lapply(tmp, function(x){as.data.table(x$p.value, keep.rownames = TRUE)})

pairwise.wilcox.test(indiv$occupancy_sum_T1, indiv$Primary.Lifestyle, p.adjust.method = "fdr")






### range change per primary lifestyle --------------
stat_res = permdiff(indiv, group="Primary.Lifestyle", rep=1000, metric="change_occupancy")
fwrite(stat_res$stats, "output/change_occupancy_primary_lifestyle_permutation_stats.csv")
primary_lifestyleplot_occ = ggplot(stat_res$data, aes(y = Primary.Lifestyle, x = change_occupancy, fill = atlas)) +
  geom_vline(data=stat_res$stats, aes(xintercept=group_median), col="grey")+
  geom_boxplot(varwidth = F, outlier.alpha = 0.2)+
  scale_fill_startrek()+
  scale_x_continuous(trans="log1p")+
  labs(y = "", x = "\U0394 occupancy") +
  theme(legend.position = "none") +
  facet_wrap(~atlas, ncol=4, scales="free_x") +
  geom_text(data = stat_res$stats,
            aes(x = y, y = group, label = significance),
            inherit.aes = FALSE,
            hjust = 0,
            size = 4)
ggsave(paste0(figure_folder,  "occupancy_change_per_lifestyle.png"), width=180, height = 40, unit="mm", bg="white", dpi=300)


stat_res = permdiff(indiv, group="Order1", rep=1000, metric="change_occupancy")
fwrite(stat_res$stats, "output/change_occupancy_order_permutation_stats.csv")
taxplot_occ = ggplot(stat_res$data, aes(y = Order1, x = change_occupancy, fill = atlas)) +
  geom_vline(data=stat_res$stats, aes(xintercept=group_median), col="grey")+
  geom_boxplot(varwidth = F, outlier.alpha = 0.2)+
  scale_fill_startrek()+
  scale_x_continuous(trans="log1p")+
  labs(y = "", x = "\U0394 occupancy") +
  theme(legend.position = "none") +
  facet_wrap(~atlas, ncol=4, scales="free_x") +
  geom_text(data = stat_res$stats,
            aes(x = y, y = group, label = significance),
            inherit.aes = FALSE,
            hjust = 0,
            size = 4)
#ggsave(paste0(figure_folder,  "occupancy_change_per_order.png"), width=180, height = 40, unit="mm", bg="white", dpi=300)


stat_res = permdiff(indiv, group="Habitat", rep=1000, metric="change_occupancy")
fwrite(stat_res$stats, "output/change_occupancy_habitat_permutation_stats.csv")
habplot_occ = ggplot(stat_res$data, aes(y = Habitat, x = change_occupancy, fill = atlas)) +
  geom_vline(data=stat_res$stats, aes(xintercept=group_median), col="grey")+
  geom_boxplot(varwidth = F, outlier.alpha = 0.2)+
  scale_fill_startrek()+
  scale_x_continuous(trans="log1p")+
  labs(y = "", x = "\U0394 occupancy") +
  theme(legend.position = "none") +
  facet_wrap(~atlas, ncol=4, scales="free_x") +
  geom_text(data = stat_res$stats,
            aes(x = y, y = group, label = significance),
            inherit.aes = FALSE,
            hjust = 0,
            size = 4)
#ggsave(paste0(figure_folder,  "occupancy_change_per_habitat.png"), width=180, height = 40, unit="mm", bg="white", dpi=300)


stat_res = permdiff(indiv, group="Trophic.Niche", rep=1000, metric="change_occupancy")
fwrite(stat_res$stats, "output/change_occupancy_trophic_niche_permutation_stats.csv")
trophplot_occ = ggplot(stat_res$data, aes(y = Trophic.Niche, x = change_occupancy, fill = atlas)) +
  geom_vline(data=stat_res$stats, aes(xintercept=group_median), col="grey")+
  geom_boxplot(varwidth = F, outlier.alpha = 0.2)+
  scale_fill_startrek()+
  scale_x_continuous(trans="log1p")+
  labs(y = "", x = "\U0394 occupancy") +
  theme(legend.position = "none") +
  facet_wrap(~atlas, ncol=4, scales="free_x") +
  geom_text(data = stat_res$stats,
            aes(x = y, y = group, label = significance),
            inherit.aes = FALSE,
            hjust = 0,
            size = 4)
#ggsave(paste0(figure_folder,  "occupancy_change_per_trophic_niche.png"), width=180, height = 40, unit="mm", bg="white", dpi=300)



plot_grid(primary_lifestyleplot_occ, taxplot_occ, habplot_occ, trophplot_occ, ncol=1, 
          labels = "auto",  
          label_fontface = "bold", label_size = 10, rel_heights = c(30, 70, 40, 40))
ggsave(paste0(figure_folder,  "figS10.png"), width=200, height = 30+65+40+45+10, unit="mm", bg="white", dpi=300)





















# same but for presence absence
# stat_res = permdiff(indiv, group="Primary.Lifestyle", rep=1000, metric="change_pres.abs")
# fwrite(stat_res$stats, "output/change_pres.abs_habitat_permutation_stats.csv")
# ggplot(stat_res$data, aes(y = Primary.Lifestyle, x = change_pres.abs, fill = atlas)) +
#   geom_vline(data=stat_res$stats, aes(xintercept=group_median), col="grey")+
#   geom_boxplot(varwidth = F, outlier.alpha = 0.2)+
#   scale_fill_startrek()+
#   scale_x_continuous(trans="log1p")+
#   labs(y = "", x = "\U0394 presence grids relative to 1^st^ sampling period") +
#   theme(legend.position = "none") +
#   facet_wrap(~atlas, ncol=4, scales="free_x") +
#   geom_text(data = stat_res$stats,
#             aes(x = y, y = group, label = significance),
#             inherit.aes = FALSE,
#             hjust = 0,
#             size = 4)
# ggsave(paste0(figure_folder,  "presence_change_per_lifestyle.png"), width=180, height = 40, unit="mm", bg="white", dpi=300)











# Balanced ranges -----

# divide species pair
chan$sp1 = gsub("\\|.*", "", chan$species_pair)
chan$sp2 = gsub(".*\\|", "", chan$species_pair)
# add ranges for species from T1, for both species (occupancy)
chan = merge(chan, unique(range_change[endYear==1, .(atlas, scientificName, mean.psi_sum)]), all.x=TRUE, 
             by.x=c("sp1", "atlas"), by.y=c("scientificName", "atlas"))
chan = merge(chan, unique(range_change[endYear==1, .(atlas, scientificName, mean.psi_sum)]), all.x=TRUE, 
             by.x=c("sp2", "atlas"), by.y=c("scientificName", "atlas"))
chan = merge(chan, unique(range_change[endYear==1, .(atlas, scientificName, pres.abs_sum)]), all.x=TRUE, 
              by.x=c("sp1", "atlas"), by.y=c("scientificName", "atlas"))
chan = merge(chan, unique(range_change[endYear==1, .(atlas, scientificName, pres.abs_sum)]), all.x=TRUE, 
             by.x=c("sp2", "atlas"), by.y=c("scientificName", "atlas"))

# mark pairs with species  < 50 grids (this is based on presence absence records in T1 to keep it straight forward)
chan[, small_ranges:= pres.abs_sum.x<50 | pres.abs_sum.y<50, ]
# mark pairs with unbalanced range ratios
chan[,range_ratio_mean.psi:=min(c(mean.psi_sum.x, mean.psi_sum.y))/max(c(mean.psi_sum.x,mean.psi_sum.y)), by=.(species_pair, atlas)]
# the minimum number of occ cells in a pair (both species)
chan[,min_pair_range:=min(c(mean.psi_sum.x,mean.psi_sum.y)), by=.(species_pair, atlas)]


### plot 
tmp = unique(chan[, .(delta_cor_obs, range_ratio_mean.psi, atlas, min_pair_range)])
ggplot(tmp, aes(y=abs(delta_cor_obs), x=range_ratio_mean.psi, col=min_pair_range))+
  geom_point(alpha=0.05)+
  geom_smooth(data = tmp[range_ratio_mean.psi<0.25,])+
  geom_smooth(data = tmp[range_ratio_mean.psi>= 0.25 & range_ratio_mean.psi<0.5,])+
  geom_smooth(data = tmp[range_ratio_mean.psi>= 0.5 & range_ratio_mean.psi<0.75,])+
  geom_smooth(data = tmp[range_ratio_mean.psi>=0.75,])+
  labs(x="range size ratio", y="Absolute \U0394 Spearman's \U03C1")+
  scale_color_viridis_c(name = "Sum occupancy probablity\n (smaller species)")+
  facet_wrap("atlas")+
  stat_cor(data = tmp[range_ratio_mean.psi<0.25,], method = "spearman", size=2, p.accuracy = 0.001, cor.coef.name = "rho",
           label.sep = "\n", col="black",label.x = 0.05)+
  stat_cor(data = tmp[range_ratio_mean.psi>= 0.25 & range_ratio_mean.psi<0.5,], method = "spearman", size=2, p.accuracy = 0.001, cor.coef.name = "rho",
           label.sep = "\n", col="black",label.x = 0.3)+
  stat_cor(data = tmp[range_ratio_mean.psi>= 0.5 & range_ratio_mean.psi<0.75,], method = "spearman", size=2, p.accuracy = 0.001, cor.coef.name = "rho",
           label.sep = "\n", col="black",label.x = 0.55)+
  stat_cor(data = tmp[range_ratio_mean.psi>=0.75,], method = "spearman", size=2, p.accuracy = 0.001, cor.coef.name = "rho",
           label.sep = "\n", col="black",label.x = 0.8)+
  theme(legend.position = "bottom", legend.key.width = unit(12, "mm"))+
  geom_vline(xintercept = c(0.25,0.5,0.75), col="grey")
ggsave(paste0(figure_folder,  "delta_cor_obs_VS_range_size_ratio.png"), width=6.5, height = 6, bg="white")

# quantify variation:
cor.test(abs(tmp$delta_cor_obs), tmp$range_ratio_mean.psi, method="s")
tapply(tmp, tmp$atlas, function(x){cor.test(x$delta_cor_obs, x$range_ratio_mean.psi, method='s')})


# get Fig 4 for more range_ratio>=0.15 & small_ranges==FALSE (SI figure)
sub = chan[range_ratio_mean.psi>=0.2 & small_ranges==FALSE,]


# get average changes for species, only from balanced pairs
indiv_br = list()
at = sort(unique(sub$atlas))
for(a in 1:length(at)){
  tmp = sub[atlas==at[a],]
  indiv_int = data.table(species = sort(unique(c(tmp$sp1, tmp$sp2))), 
                         delta_cor_obs_median = NA,
                         delta_c_obs_median = NA,
                         atlas = at[a],
                         n_pairs = NA)
  for(i in 1:nrow(indiv_int)){
    id = grep(indiv_int$species[i], tmp$species_pair)
    indiv_int$delta_cor_obs_median[i] = median(tmp$delta_cor_obs[id], na.rm=TRUE)
    indiv_int$delta_cor_obs_sum[i] = sum(tmp$delta_cor_obs[id], na.rm=TRUE)
    indiv_int$delta_cor_obs_mean[i] = mean(tmp$delta_cor_obs[id], na.rm=TRUE)
    indiv_int$delta_cor_obs_sd[i] = sd(tmp$delta_cor_obs[id], na.rm=TRUE)
    indiv_int$cor_obs_T1_median[i] = median(tmp$cor_obs_T1[id], na.rm=TRUE)
    indiv_int$delta_c_obs_median[i] = median(tmp$delta_c_obs[id], na.rm=TRUE)
    indiv_int$n_pairs[i] = length(id)
  }
  indiv_br[[a]] = indiv_int
  cat(a, i,  "\r")
}
indiv_br = rbindlist(indiv_br)

# add avonet
indiv_br = merge(indiv_br, 
                 avonet[,.SD, .SDcols = !c('Mass.Refs.Other', 'Reference.species', 'Traits.inferred')],
                 by="species", all.x=TRUE)

# check changes for Primary.Lifestyle groups
stat_res = permdiff(indiv_br, group="Primary.Lifestyle", rep=1000, metric="delta_cor_obs_median")
fwrite(stat_res$stats, "output/primarylifestyle_permutation_stats_nosmall.csv")

figs7 <- ggplot(stat_res$data, aes(y = Primary.Lifestyle, x = delta_cor_obs_median, fill = atlas)) +
  geom_vline(data=stat_res$stats, aes(xintercept = group_median), col = "grey") +
  geom_boxplot(varwidth = T, outlier.alpha = 0.2) +
  scale_fill_startrek() +
  labs(y = "", x = "Species median \u0394 Spearman's \u03C1 (balanced ranges)") +
  theme(legend.position = "none") +
  facet_wrap(~atlas, ncol = 4, scales="free_x") +
  geom_text(data = stat_res$stats,
            aes(x = y, y = group, label = significance),
            inherit.aes = FALSE,
            hjust = 0,
            size = 4)
ggsave(paste0(figure_folder,  "delta_cor_obs_per_Primary.Lifestyle_larger_ranges.png"), width=200, height = 40, unit="mm", bg="white", dpi=300)

plot_grid(scatter_delta_occ, range_diff_t1, figs7,
          ncol=1, labels="auto", rel_heights = c(60,40,40))
ggsave(paste0(figure_folder,  "occupancy_affects_and_balanced_ranges.png"),
       width=200, height = 62+40+40, unit="mm", bg="white", dpi=300)


rm(indiv_int, indiv_br)











# Variability in co-occurrence change within species --------------------------
## SI figure

# get ranks
indiv2[, median_delta_cor_obs:=median(delta_cor_obs), by=.(species, atlas)]
tmp = unique(indiv2[,.(species, median_delta_cor_obs, atlas)])
tmp[, rank:=rank(median_delta_cor_obs), by=.(atlas)]

indiv2 = merge(indiv2, tmp, by=c('atlas', 'species'))

# add ranges
indiv2[, min_range:=min(delta_cor_obs), by=.(species, atlas)]
indiv2[, max_range:=max(delta_cor_obs), by=.(species, atlas)]

# add interquartile ranges
indiv2[, min_IQR:=quantile(delta_cor_obs, 0.25), by=.(species, atlas)]
indiv2[, max_IQR:=quantile(delta_cor_obs, 0.75), by=.(species, atlas)]

ggplot(indiv2, aes(xmin = min_range, xmax = max_range, y=rank))+
  geom_linerange(col="grey")+
  geom_linerange(aes(xmin = min_IQR, xmax = max_IQR, col=atlas))+
  scale_color_startrek(guide="none") +
  geom_vline(xintercept=0, col="white")+
  geom_point(aes(x=median_delta_cor_obs.x), col="black", size=0.5)+
  labs(y="Species ordered by median co-occurrence change", x="\U0394 Spearman's \U03C1")+
  facet_wrap('atlas', scales="free_y")

ggsave(paste0(figure_folder,  "intraspecific_variability.png"), width=200, height = 160, unit="mm", bg="white", dpi=300)














