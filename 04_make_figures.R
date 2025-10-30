# Script: 04_make_figures.R
# Author: Melanie Tietje
# Email: tietje@fzp.czu.cz
# GitHub: @Eryops1
# Last Modified: 2025-10-29
# Purpose: Make figures
# Output: Figures, tables
# Notes: groundhog will ensure the exact same R packages will be used and create
#        a library on first run, which might take a while
#         

# Load libraries ----------------------------------------------------------

rm(list=(ls())) # clear workspace
library(groundhog)
pkgs <- c("data.table",
          "ggplot2",
          "cowplot",
          "ggsci",
          "ggtext")
groundhog.library(pkgs, "2025-10-25")

theme_set(theme_classic(base_size = 9)+
            theme(strip.background = element_blank())+
            theme(axis.title = ggtext::element_markdown())+
            theme(legend.title = ggtext::element_markdown()))



# functions --------------------------------------------------------------
sort_pairs_vectorized <- function(x) {
  # a function that sorts species pair names string into alphabetic order
  sapply(x, function(item) {
    spl <- unlist(strsplit(item, split = "|", fixed = TRUE))
    paste0(sort(spl), collapse = "|")
  })
}









# Read data ----------------------------------------------------------------


fin = readRDS("data/processed_spass.rds")
fin[, over:=ifelse(all(overlap==FALSE), FALSE, TRUE), by=species_pair]

# add atlas names
fin$atlas = as.character(fin$dataset_id)
fin$atlas = gsub("26", "Europe", fin$atlas)
fin$atlas = gsub("5", "Czechia", fin$atlas)
fin$atlas = gsub("6", "New York", fin$atlas)
fin$atlas = gsub("17", "New Zealand", fin$atlas)
table(fin$atlas)



# timeline ----------------------------------------------------------------
# number of species and pairs in which years

temp = unique(fin[,.(atlas, time_bin, over, species_pair)])
temp = temp[over==TRUE,]
temp$years = temp$time_bin
species = temp[,unique(unlist(strsplit(species_pair, split = "|", fixed = TRUE))), by=atlas]
species = species[,length(V1), by=atlas]
tmp = data.table(table(temp$atlas[temp$time_bin=="T1"]))
setnames(tmp, old=c('V1', 'N'), new=c('atlas', 'N pairs'))

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

ggplot(temp, aes(x=year, y=atlas, col=atlas, size=V1))+
  geom_point(pch=15)+
  scale_color_startrek(guide='none')+
  scale_size_continuous(guide='none')+
  labs(y="")+theme(panel.background = element_blank())
ggsave('figures/timeline_data.svg', height=1.5, width=2.4, bg = NULL)

# Flo says time span covered is ceiling(mean(years(atlas2))) - ceiling(mean(years(atlas1)))
mean(ceiling(mean(2014:2017)) - ceiling(mean(1985:1989)), # CZ
ceiling(mean(2000:2005)) - ceiling(mean(1980:1985)), # NY
ceiling(mean(1999:2004)) - ceiling(mean(1969:1979)), # NZ
ceiling(mean(2013:2017)) - ceiling(mean(1972:1995)) # EU
)







# STATS --------------

stat = fin[over==TRUE,]
stat[,status_T1:=cut(cor_ses[time_bin=="T1"], breaks=c(-Inf,-1.96,1.96,Inf), labels=c("segregated", "neutral", "aggregated")), 
by=.(atlas, species_pair, scaleID)]
# stat[,status_T3:=cut(cor_ses[time_bin=="T3"], breaks=c(-Inf,-1.96,1.96,Inf), labels=c("segregated", "neutral", "aggregated")), 
#      by=.(atlas, species_pair, scaleID)]
# stat[,status_T2:=cut(cor_ses[time_bin=="T2"], breaks=c(-Inf,-1.96,1.96,Inf), labels=c("segregated", "neutral", "aggregated")), 
#      by=.(atlas, species_pair, scaleID)]
stat[,status_Tlast:=cut(cor_ses[time_bin==tail(time_bin,1)], breaks=c(-Inf,-1.96,1.96,Inf), labels=c("segregated", "neutral", "aggregated")), 
    by=.(atlas, species_pair, scaleID)]
View(stat[,.(species_pair, atlas, cor_ses, time_bin, status_T1, status_Tlast)])
stat = unique(stat[,.(species_pair, atlas, status_T1, status_Tlast)])

# proportions of segregated etc....
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

fwrite(dt[, .(atlas, `co-occurrence_status`, sampling_period, years, N, prop, prop_per_group, prop_per_group_range, prop_change)], "output/status_stats1.csv")


ggplot(dt, aes(y=prop, x=`co-occurrence_status`, fill=`co-occurrence_status`))+
  geom_bar(stat = "identity")+
  facet_grid(sampling_period~atlas)

# these are not the same numbers as for the circlize plot, why - is that true?








## raw ----
ggplot(fin[over==TRUE,], aes(x=c_obs, fill=atlas))+
  geom_histogram(bins=80)+
  scale_y_sqrt()+
  stat_central_tendency(type="median", alpha=0.5, linetype=2)+
  facet_grid(time_bin~atlas, scale="free_y")+
  scale_fill_startrek(alpha=0.7, guide="none")+
  labs(x="*C*-score", y="")
  ggsave(paste0("figures/", "c_obs_histogram.png"), width=5.5, height = 5, bg="white")
  
# ggplot(fin[over==TRUE,], aes(x=cor_obs, fill=atlas))+
#   geom_histogram(bins=40)+
# #  scale_y_sqrt()+
#   stat_central_tendency(type="median", alpha=0.5, linetype=2)+
#   facet_grid(time_bin~atlas, scale="free_y")+
#   scale_fill_startrek(alpha=0.7)+
#   labs(x="Spearman's \U03C1", y="")+
#   theme(legend.position = "none")
#   ggsave(paste0("figures/", "cor_obs_histogram.png"), width=5.5, height = 5, bg="white")  

# FIG 1
# ggplot(fin[over==TRUE,], aes(x=cor_obs, y=atlas, fill=atlas, colour = atlas, lty=time_bin))+
# #  geom_density_ridges(quantile_lines = T, quantiles = 2, vline_width=0.5, vline_linetype=1,
# #                      alpha=0.2, scale=1)+
#   # geom_density_ridges(stat="binline", vline_width=0.5, vline_linetype=2, 
#   #                     alpha=0.2, scale=1)+
#   geom_vline(xintercept = 0, lty=1, col="grey50", lwd=0.2)+
#   scale_fill_startrek(guide="none")+
#   scale_color_startrek(guide="none")+
#   scale_linetype_discrete("sampling period")+
#   coord_cartesian(ylim=c(1.5,4.5))+
#   labs(y="", x="Spearman's \U03C1")+
#   theme(legend.position = "bottom")
#  ggsave(paste0("figures/", "cor_obs_ridges.png"), width=5, height = 4, bg="white")
  
fig1dt = fin[!(time_bin=="T2"& atlas=="Czechia"),]
fig1dt$time_bin[fig1dt$time_bin=="T3"] = "T2"
ggplot(fig1dt[over==TRUE,], aes(x=cor_obs))+
  geom_histogram(binwidth = 0.1, fill="grey80")+ # aes(y = ..density..), 
  #geom_density(alpha=0.1, aes(fill=atlas, colour = atlas), bw=0.05)+
  facet_grid(time_bin~atlas, scales="free_y")+
  geom_vline(xintercept = 0, lty=1, col="grey40", lwd=0.3)+
  stat_central_tendency(type="median", aes(col=atlas), lty=2, lwd=0.3)+
  scale_fill_startrek(guide="none")+
  scale_color_startrek(guide="none")+
  labs(y="count", x="Co-occurrence (Spearman's \U03C1)")+
  theme(legend.position = "none")
ggsave(paste0("figures/", "cor_obs_histogram.pdf"), width=6, height = 2.5, bg="white")

fig1dt$x_bin <- cut(
  fig1dt$cor_ses,
  breaks = c(-Inf, -1.96, 1.96, Inf),
  labels = c("Seg (< -1.96)", "Neutral (-1.96 to 1.96)", "Agg (> 1.96)"),
  ordered_result = TRUE   # ensures ordering
)
ggplot(fig1dt[over==TRUE,], aes(x=cor_obs))+
  geom_histogram(binwidth = 0.1, aes(fill=x_bin), bins=160)+ # aes(y = ..density..), , fill="grey80"
  facet_grid(time_bin~atlas, scales="free_y")+
  geom_vline(xintercept = 0, lty=1, col="grey40", lwd=0.3)+
  stat_central_tendency(type="median", col="black", lty=2, lwd=0.3)+
  #scale_fill_startrek(guide="none")+
  #grid.col = c("segregated"="#1398E9", "neutral"="grey", "aggregated"="#E96513")
  scale_fill_manual(name = "", values = c("#1398E9", "grey", "#E96513"), labels=c("segregated", "neutral", "aggregated"))+
  scale_color_startrek(guide="none")+
  labs(y="count", x="Co-occurrence (Spearman's \U03C1)")+
  theme(legend.position = "bottom")
ggsave(paste0("figures/", "cor_obs_histogram2.png"), width=6, height = 2.5, bg="white")



# ggplot(fin[over==TRUE & cor_obs_p<0.05,], aes(x=cor_obs, y=atlas, fill=atlas, colour = atlas, lty=time_bin))+
#   geom_density_ridges(quantile_lines = T, quantiles = 2, vline_width=0.5, vline_linetype=1,
#                       alpha=0.2, scale=1)+
#   geom_vline(xintercept = 0, lty=1, col="grey50", lwd=0.2)+
#   scale_fill_startrek(guide="none")+
#   scale_color_startrek(guide="none")+
#   scale_linetype_discrete("sampling period")+
#   coord_cartesian(ylim=c(1.5,4.5))+
#   labs(y="", x="Spearman's \U03C1")+
#   theme(legend.position = "bottom")
# ggsave(paste0("figures/", "cor_obs_ridges_significant.png"), width=4.5, height = 2.5, bg="white")
  

# ggplot(fin[over==TRUE,], aes(x=c_obs, y=atlas, fill=atlas, colour = atlas, lty=time_bin))+
#   geom_density_ridges(quantile_lines = F, quantiles = 2, vline_width=0.5, vline_linetype=1,
#                       alpha=0.2, scale=1)+
#   scale_fill_startrek(guide="none")+
#   scale_color_startrek(guide="none")+
#   scale_linetype_discrete("sampling period")+
#   coord_cartesian(ylim=c(1.5,4.5))+
#   theme(legend.position = "bottom")+
#   labs(y="", x="*C*-score")
# ggsave(paste0("figures/", "c_obs_ridges.png"), width=4.5, height = 2.5, bg="white")

ggplot(fin[over==TRUE,], aes(x=c_obs, y=cor_obs, fill=atlas, color=atlas))+
  geom_point(alpha=0.02, shape = 21)+
  stat_cor(method = "spearman", label.x = 0.5, label.y = -0.9, size=3, p.accuracy = 0.001, cor.coef.name = "rho")+
  facet_grid(time_bin~atlas)+
  geom_smooth(colour = "grey40")+
  scale_fill_startrek(guide="none")+
  scale_color_startrek(guide="none")+
  labs(y="Spearman's \U03C1", x="*C*-scores")
ggsave(paste0("figures/", "c_obs_cor_obs_scatterplot.png"), width=10, height = 6, bg="white", dpi = 150)  

### stats ----
psych::describeBy(data=fin, cor_obs ~ atlas + time_bin, mat=TRUE, digits=2)
pairwise.wilcox.test(fin$cor_obs[fin$dataset_id==5], fin$time_bin[fin$dataset_id==5])
pairwise.wilcox.test(fin$cor_obs[fin$dataset_id==6], fin$time_bin[fin$dataset_id==6])
pairwise.wilcox.test(fin$cor_obs[fin$dataset_id==17], fin$time_bin[fin$dataset_id==17])
pairwise.wilcox.test(fin$cor_obs[fin$dataset_id==26], fin$time_bin[fin$dataset_id==26])




# compare to PSD matrix (cor() values (CZ only))
psd= readRDS("PSD_matrix_atlas5.rds")
pw = fig1dt[over==TRUE & atlas=="Czechia",]
hist(as.dist(psd[[1]]))
hist(pw[dataset_id==5 & time_bin=="T1", cor_obs])

# turn dist into DT
psd_dt <- as.data.table(as.table(psd[[1]]))
setnames(psd_dt, c("sp1", "sp2", "value"))
psd_dt <- psd_dt[sp1 < sp2]
psd_dt[, species_pair := paste(sp1, sp2, sep = "|")]
psd_dt <- psd_dt[, .(species_pair, value)]

matrixcalc::is.positive.definite(psd[[1]])   # strict PD


# match order
psd_dt$species_pair = sort_pairs_vectorized(psd_dt$species_pair)
pw$species_pair = sort_pairs_vectorized(pw$species_pair)

table(psd_dt$species_pair %in% pw$species_pair) # ok
table(pw$species_pair %in% psd_dt$species_pair) # ok

tes = merge(psd_dt, pw[time_bin=="T1",.(species_pair, cor_obs, cor_obs_p)], )
ggplot(tes, aes(x=value, y=cor_obs, col=cor_obs_p<0.05))+
  geom_point(alpha=0.2)+
  stat_cor(method = "spearman", size=3, p.accuracy = 0.001, cor.coef.name = "rho")
tes[,diff:=value-cor_obs,]
hist(tes$diff, border=NA, xlab="matrix correlation - pairwise correlation")
abline(v=0, lty=2)







## Z-scores ----

# ggplot(fin[over==TRUE,], aes(x=cor_ses, fill=atlas, colour = atlas, lty=time_bin))+
#   geom_density(alpha=0.2)+
#   geom_vline(xintercept = 0, lty=1, col="grey50", lwd=0.2)+
#   scale_fill_startrek(guide="none")+
#   scale_color_startrek(guide="none")+
#   scale_linetype_discrete("sampling period")+
#   facet_grid("atlas", scales="free")+
#   labs(y="", x="*Z*-score Spearman's \U03C1")+
#   theme(legend.position = "bottom")

# ggplot(fin[over==TRUE,], aes(x=cor_ses, y=atlas, fill=atlas, colour = atlas, lty=time_bin))+
#   geom_density_ridges(quantile_lines = T, quantiles = 2, vline_width=0.5, vline_linetype=1,
#                       alpha=0.2, scale=1)+
#   geom_vline(xintercept = 0, lty=1, col="grey50", lwd=0.2)+
#   scale_fill_startrek(guide="none")+
#   scale_color_startrek(guide="none")+
#   scale_linetype_discrete("sampling period")+
#   coord_cartesian(ylim=c(1.5,4.5))+
#   labs(y="", x="*Z*-score Spearman's \U03C1")+
#   theme(legend.position = "bottom")
# ggsave(paste0("figures/", "cor_Zscore_ridges.png"), width=5, height = 4, bg="white")

# fig1dt = fin[!(time_bin=="T2"& atlas=="Czechia"),]
# fig1dt$time_bin[fig1dt$time_bin=="T3"] = "T2"
# ggplot(fig1dt[over==TRUE,], aes(x=cor_ses))+
#   geom_histogram(binwidth = 4, fill="grey80")+ # aes(y = ..density..), 
#   #geom_density(alpha=0.1, aes(fill=atlas, colour = atlas), bw=0.05)+
#   facet_grid(time_bin~atlas, scales="free_y")+
#   geom_vline(xintercept = 0, lty=1, col="grey40", lwd=0.3)+
#   stat_central_tendency(type="median", aes(col=atlas), lty=2)+
#   scale_fill_startrek(guide="none")+
#   scale_color_startrek(guide="none")+
#   labs(y="count", x="Co-occurrence (_Z_-score Spearman's \U03C1)")+
#   theme(legend.position = "none")
# ggsave(paste0("figures/", "cor_ses_histogram.png"), width=6, height = 2.5, bg="white")

fig1dt$x_bin <- cut(
  fig1dt$cor_ses,
  breaks = c(-Inf, -1.96, 1.96, Inf),
  labels = c("Seg (< -1.96)", "Neutral (-1.96 to 1.96)", "Agg (> 1.96)"),
  ordered_result = TRUE   # ensures ordering
)
zhist = ggplot(fig1dt[over==TRUE,], aes(x=cor_ses))+
  geom_histogram(aes(fill=x_bin), bins=40)+ # aes(y = ..density..), , fill="grey80"
  facet_grid(time_bin~atlas, scales="free_y")+
  geom_vline(xintercept = 0, lty=1, col="grey40", lwd=0.3)+
  stat_central_tendency(type="median", col="black", lty=2, lwd=0.3)+
  #scale_fill_startrek(guide="none")+
  #grid.col = c("segregated"="#1398E9", "neutral"="grey", "aggregated"="#E96513")
  scale_fill_manual(name = "", values = c("#1398E9", "grey", "#E96513"), labels=c("segregated", "neutral", "aggregated"))+
  scale_color_startrek(guide="none")+
  labs(y="count", x="Co-occurrence (_Z_-score Spearman's \U03C1)")+
  theme(legend.position = "bottom")
ggsave(paste0("figures/", "cor_ses_histogram2.png"), width=6, height = 2.5, bg="white")


# test = fin[!is.infinite(c_ses),]
# ggplot(test[over==TRUE,], aes(x=c_ses, y=atlas, fill=atlas, colour = atlas, lty=time_bin))+
#   geom_density_ridges(quantile_lines = F, quantiles = 2, vline_width=0.5, vline_linetype=2, scale=0.9, 
#                       alpha=0.2)+
#   scale_fill_manual(values=startrek, guide="none")+
#   scale_color_manual(values=startrek, guide="none")+
#   scale_linetype_discrete("sampling period")+
#   coord_cartesian(ylim=c(1.5,4.2))+
#   labs(y="", x="*Z*-score *C*-score")+
#   theme(legend.position = "bottom")
# ggsave(paste0("figures/", "c_Zscore_ridges.png"), width=5, height = 4, bg="white")

ggplot(fin[over==TRUE,], aes(x=c_ses, y=cor_ses, col=atlas))+
  geom_point(alpha=0.02)+
  stat_cor(method = "spearman", label.x = 0.5, label.y = -0.9, size=3, p.accuracy = 0.001, cor.coef.name = "rho")+
  facet_grid(time_bin~atlas, scales="free")+
  geom_smooth(color="grey40")+
  scale_color_startrek(guide=FALSE)+
  labs(y="*Z*-score Spearman's \U03C1", x="*Z*-score *C*-score")
ggsave(paste0("figures/", "c_Zscore_cor_Zscore_scatterplot.png"), width=8, height = 7, bg="white")  


### stats -----
psych::describeBy(data=fin, cor_ses ~ atlas + time_bin, mat=TRUE, digits=2)








## raw vs Z-scores ----
corscat = ggplot(fig1dt, aes(x=cor_obs, y=cor_ses, col=atlas))+
  geom_point(alpha=0.01)+
  stat_cor(method = "spearman", label.x = 0.2, label.y = -25, size=3, 
           label.sep = "\n", p.accuracy = 0.001, cor.coef.name = "rho", col="black")+
  facet_grid(time_bin~atlas)+
  scale_color_startrek(guide=FALSE)+
  labs(y="*Z*-score Spearman's \U03C1", x="Spearman's \U03C1")
ggsave(paste0("figures/", "cor_obs_cor_Zscore_scatterplot.png"), width=6, height = 4.5, bg="white", dpi = 300)  

ggplot(fin[over==TRUE,], aes(x=c_obs, y=c_ses, col=atlas))+
  geom_point(alpha=0.01)+
  stat_cor(method = "spearman", label.x = 0.6, label.y = -30, size=3, p.accuracy = 0.001, cor.coef.name = "rho")+
  facet_grid(time_bin~atlas)+
  scale_color_startrek(guide=FALSE)+
  labs(y="*Z*-score *C*-score", x="*C*-score")
ggsave(paste0("figures/", "c_obs_c_Zscore_scatterplot.png"), width=9, height = 7, bg="white", dpi = 150)  

plot_grid(zhist, corscat, ncol=1, labels = c("a)", "b)"), label_fontface = "plain", label_size = 10)
ggsave(paste0("figures/", "figS1.png"), width=6, height = 6, bg="white")

# # all datasets,  scale 1
# ggplot(res[scaleID=="S1" | 
#              (scaleID=="S4" & dataset_id==36) | 
#              (scaleID=="S2" & dataset_id==6),], aes(x=c_ses))+
#   geom_histogram(bins=80)+
#   geom_vline(xintercept = 0, lty=1, col="grey", lwd=0.1)+
#   geom_vline(aes(xintercept = median), lty=1, col="red", lwd=0.1)+
#   xlab("Z-score of C-scores")+
#   scale_y_sqrt()+
#   facet_grid(dataset_id~time_bin, scale="free_y")+
#   coord_cartesian(expand = F)














# DYNAMICS  ---------------------------------------------------------

chan = unique(fin[,.(over, species_pair, atlas, scaleID, 
                     delta_cor_obs, delta_c_obs, delta_c_obs_z, delta_cor_obs_z)])
chan = chan[over==TRUE,]
tapply(chan$species_pair, chan$atlas, function(x){length(unique(x))})



## raw ----

# ggplot(chan, aes(x=delta_cor_obs, y=atlas, col=atlas, fill=atlas))+
#   geom_density_ridges(quantile_lines = T, quantiles = 2, vline_width=0, vline_linetype=2)+
#   geom_vline(xintercept = 0, lty=1, col="grey20", lwd=0.1)+
#   labs(y="", x="\U0394 Spearman's \U03C1")+
#   stat_central_tendency(type="median", aes(col=atlas), lty=2)+
#   scale_color_startrek(guide='none', alpha=0.5)+
#   scale_fill_startrek(guide='none', alpha=0.5)+
#   coord_cartesian(ylim=c(1.5,5.2))
# ggsave(paste0("figures/", "delta_cor_obs_ridges.png"), width=4.5, height = 2.5, bg="white")

# V2
chan[,delta_cor_obs_median:=median(delta_cor_obs, na.rm=TRUE), by=atlas]
chan[,delta_cor_obs_mean:=mean(delta_cor_obs, na.rm=TRUE), by=atlas]
statdat = unique(chan[,.(delta_cor_obs_median, delta_cor_obs_mean,atlas)])
statdat$delta_cor_obs_median = round(statdat$delta_cor_obs_median,3)
statdat$delta_cor_obs_mean = round(statdat$delta_cor_obs_mean,3)
statdat$label_med = paste0('median:\n', statdat$delta_cor_obs_median)
statdat$label_m = paste0('mean:\n', statdat$delta_cor_obs_mean)
(fig2 = ggplot(chan, aes(x=delta_cor_obs, col=atlas, fill=atlas))+
  geom_density()+
  facet_wrap(~atlas, ncol=4)+
  #scale_y_sqrt()+
  labs(y="density", x="Co-occurrence change\n(\U0394 Spearman's \U03C1)")+
  stat_central_tendency(type="mean", aes(col=atlas), lty=2)+
  scale_color_startrek(guide='none')+
  scale_fill_startrek(guide='none', alpha=0.5)+
  geom_text(data=statdat, aes(label=label_m), y=2, x=0.5, col="black", size=2.5))
ggsave(paste0("figures/", "delta_cor_obs_density.png"), width=6, height = 1.7, bg="white", dpi=300)

# # V2 poster
# ggplot(chan, aes(x=delta_cor_obs, col=atlas, fill=atlas))+
#   geom_density()+
#   facet_wrap(~atlas, ncol=2)+
#   #scale_y_sqrt()+
#   labs(y="density", x="Co-occurrence change\n(\U0394 Spearman's \U03C1)")+
#   stat_central_tendency(type="mean", aes(col=atlas), lty=2)+
#   scale_color_startrek(guide='none')+
#   scale_fill_startrek(guide='none', alpha=0.5)+
#   geom_text(data=statdat, aes(label=label_m), y=2, x=0.5, col="black", size=2.5)
# ggsave(paste0("figures/", "delta_cor_obs_density_2x2.png"), width=3, height = 3, bg="white", dpi=600)


# C score
chan[,delta_c_obs_median:=median(delta_c_obs, na.rm=TRUE), by=atlas]
chan[,delta_c_obs_mean:=mean(delta_c_obs, na.rm=TRUE), by=atlas]
statdat = unique(chan[,.(delta_c_obs_median, delta_c_obs_mean,atlas)])
statdat$delta_c_obs_median = round(statdat$delta_c_obs_median,3)
statdat$delta_c_obs_mean = round(statdat$delta_c_obs_mean,3)
statdat$label_med = paste0('median:\n', statdat$delta_c_obs_median)
statdat$label_m = paste0('mean:\n', statdat$delta_c_obs_mean)
ggplot(chan, aes(x=delta_c_obs, col=atlas, fill=atlas))+
  geom_density()+
  facet_wrap(~atlas, ncol=4)+
  labs(y="density", x="Co-occurrence change\n(\U0394 *C*-score)")+
  stat_central_tendency(type="mean", aes(col=atlas), lty=2)+
  scale_color_startrek(guide='none')+
  scale_fill_startrek(guide='none', alpha=0.5)+
  geom_text(data=statdat, aes(label=label_m), y=2, x=0.5, col="black", size=2.5)
ggsave(paste0("figures/", "delta_c_obs_density.png"), width=6, height = 1.7, bg="white", dpi=300)

ggplot(chan, aes(x=delta_c_obs, y=delta_cor_obs))+
  geom_point(alpha=0.05)+
  stat_cor(method = "spearman", label.x = 0.1, label.y = 1.4)+
  geom_smooth()

### stats -----
psych::describeBy(data=fin, delta_cor_obs ~ atlas, mat=TRUE, digits=2, IQR=TRUE, quant=c(.25,.75))





## Z-score ----

# ggplot(chan, aes(x=delta_cor_obs_z, y=atlas, col=atlas, fill=atlas))+
#   geom_density_ridges(quantile_lines = T, quantiles = 2, vline_width=0, vline_linetype=2)+
#   geom_vline(xintercept = 0, lty=1, col="grey20", lwd=0.1)+
#   scale_color_startrek(guide=FALSE, alpha=0.5)+
#   scale_fill_startrek(guide=FALSE, alpha=0.5)+
#   coord_cartesian(ylim=c(1.5,5))+
#   labs(y="", x="*Z*-score \U0394 Spearman's \U03C1")
# ggsave(paste0("figures/", "delta_c_obs_Zscore_riges.png"), width=4.5, height = 2.5, bg="white")
# 
# ggplot(chan, aes(x=delta_cor_obs_z, fill=atlas))+
#   geom_histogram(bins=80)+
#   geom_vline(xintercept = 0, lty=1, col="grey", lwd=0.1)+
#   stat_central_tendency(type="median", alpha=0.5, lty=2)+
#   facet_grid(~atlas, scale="free_y")+
#   scale_fill_startrek(guide=FALSE, alpha=0.5)+
#   labs(y="", x="*Z*-score \U0394 Spearman's \U03C1")
# ggsave(paste0("figures/", "delta_c_obs_Zscore_histogram.png"), width=6.5, height = 2.5, bg="white")

chan[,delta_cor_obs_z_median:=median(delta_cor_obs_z, na.rm=TRUE), by=atlas]
chan[,delta_cor_obs_z_mean:=mean(delta_cor_obs_z, na.rm=TRUE), by=atlas]
statdat = unique(chan[,.(delta_cor_obs_z_median, delta_cor_obs_z_mean,atlas)])
statdat$delta_cor_obs_z_median = round(statdat$delta_cor_obs_z_median,3)
statdat$delta_cor_obs_z_mean = round(statdat$delta_cor_obs_z_mean,3)
statdat$label_med = paste0('median:\n', statdat$delta_cor_obs_z_median)
statdat$label_m = paste0('mean:\n', statdat$delta_cor_obs_z_mean)
(fig2s = ggplot(chan, aes(x=delta_cor_obs_z, col=atlas, fill=atlas))+
    geom_density()+
    facet_wrap(~atlas, ncol=4)+
    #scale_y_sqrt()+
    labs(y="density", x="Co-occurrence change (*Z*-score \U0394 Spearman's \U03C1)")+
    stat_central_tendency(type="mean", aes(col=atlas), lty=2)+
    scale_color_startrek(guide='none')+
    scale_fill_startrek(guide='none', alpha=0.5)+
    geom_text(data=statdat, aes(label=label_m), y=0.04, x=20, col="black", size=2.5))
ggsave(paste0("figures/", "delta_cor_obs_z_density.png"), width=6, height = 1.7, bg="white", dpi=300)




### stats ----
psych::describeBy(data=chan, delta_c_obs ~ atlas, mat=TRUE, digits=2)











# TRANSITIONS ---------------------------------------------------------------
library(circlize)

## circlize --------------
tmp = fin[over==TRUE,]
tmp[,status_T1:=cut(cor_ses[time_bin=="T1"], breaks=c(-Inf,-1.96,1.96,Inf), labels=c("segregated", "neutral", "aggregated")), 
    by=.(atlas, species_pair, scaleID)]
tmp[,status_Tlast:=cut(cor_ses[time_bin==tail(time_bin,1)], breaks=c(-Inf,-1.96,1.96,Inf), labels=c("segregated", "neutral", "aggregated")), 
    by=.(atlas, species_pair, scaleID)]
# 
#View(tmp[,.(species_pair, time_bin, scaleID, cor_ses, dataset_id, status_T1, status_Tlast)])
tmp = unique(tmp[,.(species_pair, atlas,status_T1, status_Tlast)])
tmp[,transition:=paste(status_T1, status_Tlast, sep = "-->")]

# atli
atli = sort(unique(tmp$atlas))
grid.col = c("segregated"="#1398E9", "neutral"="grey", "aggregated"="#E96513")

png("figures/circlize.png", width =11, height=11.4, units = "cm", bg = "transparent", res = 300)
par(mfrow=c(2,2))
  for(a in 1:length(atli)){
    tmp_cz = as.data.table(table(tmp$transition[tmp$atlas==atli[a]]))
    tmp_cz[, from:=gsub("-->.*", "", V1)]
    tmp_cz[, to:=gsub(".*-->", "", V1)]
    setnames(tmp_cz, 'N', 'value')
    tmp_cz = tmp_cz[,.(from, to, value)]
    #circos.par(start.degree = 0)
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
    cats <- unique(tmp_cz$from)
    cat_props <- round(tmp_cz$value/sum(tmp_cz$value)*100, 1)
    names(cat_props) = tmp_cz$from
    
    # sum percentages per category (total across transitions)
    cat_totals <- sapply(cats, function(cat) {
      sum(cat_props[grepl(cat, names(cat_props))], na.rm = TRUE)
    })
    circos.track(track.index = 1, panel.fun = function(x, y) {
      sector <- get.cell.meta.data("sector.index")
      
      x_pos <- mean(get.cell.meta.data("xlim"))
      y_pos <- get.cell.meta.data("ylim")[2] + 1.2  # place above the sector
      if (sector %in% names(cat_totals)) {
        label <- paste0(sector, " (", format(cat_totals[sector], nsmall = 1), "%)")
        circos.text(x = x_pos, y = y_pos,
                    labels = label,
                    facing = "bending.inside", niceFacing = TRUE, cex = 0.9, col = "black")
      }
    },
    , bg.border = NA)
    title(main = atli[a], cex.main=1.1, font.main=1)
  }
dev.off()




library(ggalluvial)
tmp = fin[over==TRUE,]
tmp[,status:=cut(cor_ses, breaks=c(-Inf,-1.96,1.96,Inf), labels=c("segregated", "neutral", "aggregated")),
    by=.(atlas, species_pair, scaleID, time_bin)]
tmp = tmp[,.(time_bin, status, species_pair, atlas)]

# replace T3 with T2 for CZ!!!!!!
tmp = tmp[!(atlas=="Czechia" & time_bin=="T2"),]
tmp$time_bin[(tmp$time_bin =="T3" & tmp$atlas=="Czechia")] = "T2"


wide <- dcast(tmp[time_bin %in% c("T1", "T2")],
              species_pair + atlas ~ time_bin, value.var = "status")
transitions <- wide[, .N, by = .(T1, T2, atlas)]
transitions[,percent:=N/sum(N), by=atlas]

ggplot(transitions, aes(axis1 = T1, axis2 = T2, y = N)) +
  geom_alluvium(aes(fill = T1), width = 0/12, alpha=0.9) +
  geom_stratum(width = 1/17, fill="grey80", color = "grey40") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), angle = 90, size=2.5) +
  scale_x_discrete(limits = c("T1", "T2"), expand = c(.05, .05)) +
  scale_fill_manual(values = c("#1398E9", "grey", "#E96513"), guide="none")+
  labs(y = "count")+
  facet_wrap(~atlas, scales="free")+
  theme_void()

chord = transitions
chord$T1 = paste0(chord$T1, "_T1")
chord$T2 = paste0(chord$T1, "_T2")
grid.col = c("segregated_T1"="#1398E9", "neutral_T1"="grey", "aggregated_T1"="#E96513",
             "segregated_T1_T2"="#1398E9", "neutral_T1_T2"="grey", "aggregated_T1_T2"="#E96513")

circos.par(start.degree = 180)
chordDiagram(chord[atlas=="Europe",],
             directional = 1, direction.type = c("arrows"),
             link.arr.type= "big.arrow",
             #link.arr.length = ifelse(link.arr.type == "big.arrow", 0.02, 0.4),
             link.arr.length = 0.1, big.gap = 1,
             link.largest.ontop=F, grid.col = grid.col, 
             transparency=0.3,
             target.prop.height=0, 
             #annotationTrack = c("grid"),
             annotationTrackHeight = c(grid = 0.08, name=0.13))
chordDiagram(transitions[atlas=="Europe",],
             directional = 1, direction.type = c("arrows"),
             link.arr.type= "big.arrow",
             #link.arr.length = ifelse(link.arr.type == "big.arrow", 0.02, 0.4),
             link.arr.length = 0.1, big.gap = 20,
             link.largest.ontop=F, grid.col = grid.col, 
             transparency=0.3,
             target.prop.height=0, self.link = 1,
             #annotationTrack = c("grid"),
             annotationTrackHeight = c(grid = 0.08, name=0.13))







# CORRELATIONS ----------------------------------------------------------------


## RTM ------------

fin[, c_obs_T1:=c_obs[time_bin=="T1"], by=.(atlas, scaleID, species_pair)]
fin[, cor_obs_T1:=cor_obs[time_bin=="T1"], by=.(atlas, scaleID, species_pair)]
chan = unique(fin[,.(over, species_pair, atlas, scaleID, 
                     delta_cor_obs, delta_c_obs, delta_c_obs_z, delta_cor_obs_z,
                     c_obs_T1, cor_obs_T1)])
chan = chan[over==TRUE,]

cor.test(chan$delta_cor_obs, chan$cor_obs_T1, method="s")
cor.test(chan$delta_c_obs, chan$c_obs_T1, method="s")

cor.test(chan$delta_cor_obs_z, chan$cor_obs_T1, method="s")
cor.test(chan$delta_c_obs_z, chan$c_obs_T1, method="s")

## bins
#chan[,cor_obs_T1_bin:=cut(cor_obs_T1, breaks=5), by=atlas]

ggplot(chan, aes(x=cor_obs_T1, y=delta_cor_obs, col=atlas))+
  geom_point(alpha=0.01)+
  facet_wrap('atlas')+
  stat_cor(method = "spearman", label.x = -1, label.y = -1, size=3, 
           p.accuracy = 0.001, cor.coef.name = "rho", col="black")+
  geom_smooth(alpha=1)+
  geom_vline(xintercept = 0, col='grey') + geom_hline(yintercept = 0, col='grey')+
  scale_color_startrek(guide="none")+
  labs(y="\U0394 Spearman's \U03C1", x="Spearman's \U03C1 T1")
ggsave(paste0("figures/", "delta_c_obs_cor_obsT1_scatterplot.png"), width=6, height = 6, bg="white")

# starting correlation plays a role, but that is also an RTM since we
# can only go down from a high value

quantile(chan$delta_cor_obs[chan$atlas=="Czechia"], 0.85)
summary(chan$delta_cor_obs[chan$atlas=="Czechia"])


tmp = chan[chan$atlas=="Czechia",]
tmp_sub = tmp[which(tmp$delta_cor_obs>=quantile(tmp$delta_cor_obs, 0.95)),]
median(tmp$delta_cor_obs, na.rm=TRUE)
sd(tmp$delta_cor_obs, na.rm=TRUE)
median(tmp_sub$delta_cor_obs, na.rm=TRUE)

# simulate changes in correlation by adding same process to each start point, then cut them to whats max possible.
d = unique(data.table(species_pair = chan$species_pair, 
               cor_obs_T1 = chan$cor_obs_T1,
               delta_cor_obs = chan$delta_cor_obs))
average_change = mean(d$delta_cor_obs)
change_dist = rnorm(1000, mean=0, sd=sd(d$delta_cor_obs))
d$delta_cor_obs_scramble = sample(d$delta_cor_obs)
d$cor_obs_T2 = d$cor_obs_T1 + d$delta_cor_obs_scramble
hist(d$cor_obs_T1)
hist(d$cor_obs_T2, xlim=c(-1,1))
d$cor_obs_T2[]


# simpsons?
chan$cor_obs_T1_cut = cut(chan$cor_obs_T1, breaks=seq(-1,1,0.1))

ggplot(chan, aes(x=cor_obs_T1, y=delta_cor_obs, col=atlas, group=cor_obs_T1_cut))+
  geom_point(alpha=0.01)+
  facet_wrap('atlas')+
#  stat_cor(method = "spearman", label.x = -1, label.y = -1, size=3, 
#           p.accuracy = 0.001, cor.coef.name = "rho", col="black")+
  geom_smooth(alpha=1, method="lm")+
  geom_vline(xintercept = 0, col='grey') + geom_hline(yintercept = 0, col='grey')+
  scale_color_startrek(guide="none")+
  labs(y="\U0394 Spearman's \U03C1", x="Spearman's \U03C1 T1")


ggplot(chan, aes(y=cor_obs_T1_cut, x=delta_cor_obs, fill=atlas))+
  #geom_point(alpha=0.01)+
  #geom_histogram(alpha=0.01)+
  geom_density_ridges(quantile_lines = T, quantiles = 2, vline_width=1, vline_linetype=2)+
  #stat_binline(bins=20)+
  facet_wrap('atlas')+
  #stat_cor(method = "spearman", label.x = -1, label.y = -1, size=3, 
  #         p.accuracy = 0.001, cor.coef.name = "rho", col="black")+
  #geom_smooth(alpha=1)+
  geom_vline(xintercept = 0, col='grey40', lty=2) + #geom_hline(yintercept = 0, col='grey')+
  scale_fill_startrek(guide="none")+
  coord_cartesian(xlim=c(-1,1))+
  labs(x="\U0394 Spearman's \U03C1", y="Spearman's \U03C1 T1")

ggplot(chan, aes(x=delta_cor_obs, col=atlas))+
  geom_histogram(alpha=0.01)+
  facet_grid(cor_obs_T1_cut~atlas)+
  #stat_cor(method = "spearman", label.x = -1, label.y = -1, size=3, 
  #         p.accuracy = 0.001, cor.coef.name = "rho", col="black")+
  #geom_smooth(alpha=1)+
  geom_vline(xintercept = 0, col='grey40', lty=2) + #geom_hline(yintercept = 0, col='grey')+
  scale_color_startrek(guide="none")+
  labs(x="\U0394 Spearman's \U03C1", y="Spearman's \U03C1 T1")

ggplot(chan[atlas=="Czechia",], aes(x=delta_cor_obs))+
  geom_histogram(alpha=0.01)+
  facet_grid("cor_obs_T1_cut", scales="free")+
  coord_cartesian(xlim=c(-1,1))+
  #stat_cor(method = "spearman", label.x = -1, label.y = -1, size=3, 
  #         p.accuracy = 0.001, cor.coef.name = "rho", col="black")+
  #geom_smooth(alpha=1)+
  geom_vline(xintercept = 0, col='grey40', lty=2) + #geom_hline(yintercept = 0, col='grey')+
  scale_color_startrek(guide="none")+
  labs(x="\U0394 Spearman's \U03C1", y="Spearman's \U03C1 T1")

  
  
  

## Range size  ------------

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

# remove East Europe / Russia
# if(26 %in% dataset_id){
#   files = dir("data", pattern = "gadm41", full.names = T)
#   rem = lapply(files, vect)
#   rem = vect(rem)
#   
#   crs(rem, proj=T) == crs(grid, proj=T)
#   #new = terra::crop(grid, rem)
#   tmp = terra::intersect(grid,rem)
#   tmp$remove = 1
#   grid_new = merge(grid, tmp, all.x=T)
#   grid = grid_new[is.na(grid_new$remove),]
# }
# set scale to 2
if(6 %in% dataset_id){
  dat$scalingID[dat$datasetID==6]=2
}

# subset to species that occur in all time periods:
## step 1: remove all NAs
dat = na.omit(dat)
## step2: count years with data
dat[, alltime:=length(unique(endYear)), by=.(datasetID, scientificName)]
table(
  dat$alltime, dat$datasetID)

# delete species that do not meet required number of sample periods
cond = as.data.table(table(gsub(".*occupancy_|_.*\\.rds", "", files)))
dat = dat[datasetID==cond$V1[1] & alltime==cond$N[1] | 
            datasetID==cond$V1[2] & alltime==cond$N[2] | 
            datasetID==cond$V1[3] & alltime==cond$N[3] | 
            datasetID==cond$V1[4] & alltime==cond$N[4]]
table(dat$alltime, dat$datasetID)
any(is.na(dat$mean.psi))

# get change in occupancy probability
# change = sum occupancy last time - sum occupancy first time (?)
dat[, mean.psi_sum:=sum(mean.psi), by=.(datasetID, scientificName, endYear)]
dat[, change_occupancy:=(sum(mean.psi[endYear==max(endYear)]) - sum(mean.psi[endYear==1])) /sum(mean.psi[endYear==1]), 
    by=.(datasetID, scientificName)]

# get change in present grid cells
dat[, pres.abs_sum:=sum(pres.abs), by=.(datasetID, scientificName, endYear)]
dat[, change_pres.abs:=(sum(pres.abs[endYear==max(endYear)]) - sum(pres.abs[endYear==1])) /sum(pres.abs[endYear==1]), 
    by=.(datasetID, scientificName)]

# merge into the chan data object
dat$atlas = dat$datasetID
dat$atlas = gsub("26", "Europe", dat$atlas)
dat$atlas = gsub("5", "Czechia", dat$atlas)
dat$atlas = gsub("6", "New York", dat$atlas)
dat$atlas = gsub("17", "New Zealand", dat$atlas)

saveRDS(dat, "data/dat.rds")


# change ################
range = unique(dat[,.(endYear, mean.psi_sum, pres.abs_sum, change_occupancy, change_pres.abs, scientificName, atlas)])
avonet = fread("../AVONet/16586228/AVONET1_BirdLife.csv")
avonet$scientificName = gsub(" ", "_", avonet$Species1)
range = merge(range, 
              avonet[,.SD, .SDcols = !c('Mass.Refs.Other', 'Reference.species', 'Traits.inferred')],
              by="scientificName", all.x=TRUE)
# HABITAT
## change occupancy sum
ggplot(range, aes(x=change_occupancy, y=Habitat))+
  geom_boxplot()+
  geom_vline(xintercept = 0, lty=1, col="grey20")+
  scale_x_continuous()+
  facet_wrap("atlas")+
  labs(x="")
  theme_linedraw()

## change presences
ggplot(range, aes(x=change_pres.abs, y=Habitat))+
  geom_boxplot()+
  geom_vline(xintercept = 0, lty=1, col="grey20")+
  scale_x_continuous(transform = 'log1p')+
  facet_wrap("atlas")+
  theme_linedraw()

# LIFESTYLE
## change occupancy sum
ggplot(range, aes(x=change_occupancy, y=Primary.Lifestyle))+
  geom_boxplot()+
  geom_vline(xintercept = 0, lty=1, col="grey20")+
  scale_x_continuous()+
  facet_wrap("atlas")+
  theme_linedraw()

## change presences
ggplot(range, aes(x=change_pres.abs, y=Habitat))+
  geom_boxplot()+
  geom_vline(xintercept = 0, lty=1, col="grey20")+
  scale_x_continuous(transform = 'log1p')+
  facet_wrap("atlas")+
  theme_linedraw()


###############
tmp = unique(dat[endYear==1, .(scientificName, pres.abs_sum, atlas)])

chan$sp1 = gsub("\\|.*", "", chan$species_pair)
chan$sp2 = gsub(".*\\|", "", chan$species_pair)
chan = merge(chan, tmp, all.x=TRUE, by.x=c("sp1", "atlas"), by.y=c("scientificName", "atlas"))
chan = merge(chan, tmp, all.x=TRUE, by.x=c("sp2", "atlas"), by.y=c("scientificName", "atlas"))
chan[,range_ratio_pres.abs:=min(c(pres.abs_sum.x,pres.abs_sum.y))/max(c(pres.abs_sum.x,pres.abs_sum.y)), by=.(species_pair, atlas)]

# the minimum number of occ cells in a pair (both species)
chan[,min_pair_range:=min(c(pres.abs_sum.x,pres.abs_sum.y)), by=.(species_pair, atlas)]


### plot 
ggplot(unique(chan[, .(delta_cor_obs, range_ratio_pres.abs, atlas, min_pair_range)]), 
       aes(y=delta_cor_obs, x=range_ratio_pres.abs, col=min_pair_range))+
  geom_point(alpha=0.05)+
  geom_smooth()+
  labs(x="range size ratio", y="\U0394 Spearman's \U03C1")+
  scale_color_viridis_c(name = "N grid cells\nsmall species")+
  facet_wrap("atlas")+
  theme(legend.position = "bottom")
ggsave(paste0("figures/", "delta_cor_obs_VS_range_size_ratio.png"), width=6.5, height = 6, bg="white")
  
# ggplot(unique(chan[, .(delta_cor_obs, range_ratio_pres.abs, atlas, min_pair_range)]), 
#        aes(y=delta_cor_obs, x=min_pair_range, col=atlas))+
#   geom_point(alpha=0.05)+
#   geom_smooth(col="grey40")+
#   labs(x="N grids smaller species", y="\U0394 Spearman's \U03C1")+
#   scale_color_startrek(guide=FALSE)+
#   facet_wrap("atlas", scales="free")
# ggsave(paste0("figures/", "delta_cor_obs_VS_min_range.png"), width=6.5, height = 3, bg="white")
#   
# 
# ggplot(unique(chan[, .(delta_cor_obs, range_ratio_pres.abs, atlas, min_pair_range)]), 
#        aes(y=delta_cor_obs, x=range_ratio_pres.abs, col=atlas))+
#   geom_point(alpha=0.05)+
#   geom_smooth(col="grey40")+
#   labs(x="range size ratio", y="\U0394 Spearman's \U03C1")+
#   scale_color_startrek(guide=FALSE)+
#   facet_wrap("atlas", scales="free")
# ggsave(paste0("figures/", "delta_cor_obs_VS_range_ratio_simple.png"), width=6.5, height = 3, bg="white")
  
tmp = unique(chan[, .(delta_cor_obs, range_ratio_pres.abs, atlas, min_pair_range)])
tapply(tmp, tmp$atlas, function(x){cor.test(x$delta_cor_obs, x$range_ratio_pres.abs, method="s")})  
#tapply(tmp, tmp$atlas, function(x){cor.test(abs(x$delta_cor_obs), x$range_ratio_pres.abs, method="s")})  
  
  
  
### Larger ranges / ratios only  -------------------

## pairs with min 50 grids and max 90% of total grids
tapply(chan$pres.abs_sum.x, chan$atlas, range)
  
chan[ , balanced_range_ratio:=range_ratio_pres.abs>0.25]
chan[ , upper_range:=max(c(pres.abs_sum.x, pres.abs_sum.y)), by=atlas]
chan[, average_ranges:= pres.abs_sum.x>50 & pres.abs_sum.y>50 &
       pres.abs_sum.x<upper_range*0.9 & pres.abs_sum.y<upper_range*0.9, ]

(fig2_balance = ggplot(chan[average_ranges==TRUE & balanced_range_ratio==TRUE, ], aes(x=delta_cor_obs, y=atlas, col=atlas, fill=atlas))+
  geom_density_ridges(quantile_lines = T, quantiles = 2, vline_width=0, vline_linetype=2)+
  geom_vline(xintercept = 0, lty=1, col="grey20", lwd=0.1)+
  labs(y="", x="\U0394 Spearman's \U03C1")+
  scale_color_startrek(guide=FALSE, alpha=0.5)+
  scale_fill_startrek(guide=FALSE, alpha=0.5)+
  coord_cartesian(ylim=c(1.5,5)))
ggsave(paste0("figures/", "delta_cor_obs_ridges_large_ranges.png"), width=4.5, height = 2.5, bg="white")




saveRDS(chan, "data/chan.rds")








# Distance matrix correlation ------------------------------------------------

#source("diss_models.R")
tmp = readRDS("data/mantel_results_new.rds")
mantel_PD_res = tmp$mantel_PD_res
mantel_FD_res = tmp$mantel_FD_res
  
dd = readRDS("data/distance_matrix.rds")
pd = readRDS("data/PD_distance_matrix.rds")
fd = readRDS("data/FD_distance_matrix.rds")

## plots -----------------------

dist = data.table(correlation_change = unlist(dd),
           atlas = rep(names(dd), lengths(dd)),
           PD_distance = unlist(pd),
           FD_distance = unlist(fd))
dist$atlas = gsub("26", "Europe", dist$atlas)
dist$atlas = gsub("5", "Czechia", dist$atlas)
dist$atlas = gsub("6", "New York", dist$atlas)
dist$atlas = gsub("17", "New Zealand", dist$atlas)

ggplot(dist, aes(x=correlation_change, y=PD_distance, col=atlas))+
  geom_point(alpha=0.01, col="grey20")+
  facet_wrap(~atlas)+
  geom_smooth()+
  labs(x="\U0394 Spearman's \U03C1", y="Phylogenetic distance")+
  scale_color_startrek(guide="none")

### PD ----
perms <- as.data.table(sapply(mantel_PD_res, "[[", "perm"))
perms = melt(perms, measure.vars = names(perms), value.name = "PD_perm", variable.name = "datasetID")
stat <- data.table(statistic = sapply(mantel_PD_res, "[[", "statistic"),
                   p_value = sapply(mantel_PD_res, "[[", "signif"),
                      datasetID = names(sapply(mantel_PD_res, "[[", "statistic")))
perms <- merge(perms, stat)
perms$atlas = perms$datasetID
perms$atlas = gsub("26", "Europe", perms$atlas)
perms$atlas = gsub("5", "Czechia", perms$atlas)
perms$atlas = gsub("6", "New York", perms$atlas)
perms$atlas = gsub("17", "New Zealand", perms$atlas)
perms$metric = "PD"
# (pd_plot= ggplot(perms, aes(x=PD_perm, fill=atlas))+
#   geom_histogram(alpha=0.7)+
#   geom_vline(aes(xintercept = statistic), lty=2, col="red")+
#   scale_fill_startrek(guide="none")+
#   labs(x="Mantel test score (\U03C1)")+
#   facet_grid(~atlas, scales="free"))
# #ggsave(paste0("figures/", "mantel_PD.png"), width=6, height = 4, bg="white")
pd_tab = unique(perms[, .(datasetID, statistic, p_value, atlas)])

keep = perms

### FD -----
perms <- as.data.table(sapply(mantel_FD_res, "[[", "perm"))
perms = melt(perms, measure.vars = names(perms), value.name = "FD_perm", variable.name = "datasetID")
stat <- data.table(statistic = sapply(mantel_FD_res, "[[", "statistic"),
                   p_value = sapply(mantel_PD_res, "[[", "signif"),
                   datasetID = names(sapply(mantel_FD_res, "[[", "statistic")))
perms <- merge(perms, stat)
perms$atlas = perms$datasetID
perms$atlas = gsub("26", "Europe", perms$atlas)
perms$atlas = gsub("5", "Czechia", perms$atlas)
perms$atlas = gsub("6", "New York", perms$atlas)
perms$atlas = gsub("17", "New Zealand", perms$atlas)
perms$metric = "FD"
keep = rbind(keep, perms, use.names=FALSE)

fd_tab = unique(perms[, .(datasetID, statistic, p_value, atlas)])



# table
mantel_results_pd <- data.table(
  atlas_ID = names(mantel_PD_res),
  Matrix1 = rep("spatial association", 4),
  Matrix2 = rep("phylogenetic distance", 4),
  Mantel_r = round(sapply(mantel_PD_res, function(x) x$statistic),3),
  p_value = sapply(mantel_PD_res, function(x) x$signif),
  Permutations = sapply(mantel_PD_res, function(x) x$permutations)
)
mantel_results_fd <- data.table(
  atlas_ID = names(mantel_FD_res),
  Matrix1 = rep("spatial association", 4),
  Matrix2 = rep("functional distance", 4),
  Mantel_r = round(sapply(mantel_FD_res, function(x) x$statistic), 3),
  p_value = sapply(mantel_FD_res, function(x) x$signif),
  Permutations = sapply(mantel_FD_res, function(x) x$permutations)
)
tab = rbind(mantel_results_pd, mantel_results_fd)
write.csv(tab, "output/mantel_results.csv")



tab$atlas = as.character(tab$atlas_ID)
tab$atlas = gsub("26", "Europe", tab$atlas)
tab$atlas = gsub("5", "Czechia", tab$atlas)
tab$atlas = gsub("6", "New York", tab$atlas)
tab$atlas = gsub("17", "New Zealand", tab$atlas)
tab$significance <- ifelse(tab$p_value < 0.001, "***",
                                  ifelse(tab$p_value < 0.01, "**",
                                         ifelse(tab$p_value < 0.05, "*", "")))
tab$y_pos = rep(c(100,140,140,80),2)
tab$x_pos = tab$Mantel_r
tab$metric <- ifelse(tab$Matrix2 == "phylogenetic distance", "PD", "FD")


ggplot(keep, aes(x=PD_perm, fill=atlas))+
  geom_histogram(alpha=0.7)+
  geom_vline(aes(xintercept = statistic), lty=2, col="black")+
  scale_fill_startrek(guide="none")+
  labs(x="Mantel test score (\U03C1)", y="counts")+
  facet_grid(atlas~metric, scales="free_y")+
  geom_text(data=tab, aes(x=Mantel_r+0.01, y=y_pos, label = significance),
            inherit.aes = FALSE,
            hjust = 0,
            size = 3)

ggsave(paste0("figures/", "mantel_histograms_new.png"), width=120, height = 80, unit="mm", bg="white", dpi=300)
#ggsave(paste0("figures/", "mantel_FD.png"), width=6, height = 4, bg="white")







# SPECIES ------------------------------------------

chan = readRDS("data/chan.rds")

indiv = list()
at = sort(unique(chan$atlas))
indiv2 = list()
for(a in 1:length(at)){
  tmp = chan[atlas==at[a],]
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





# ggplot(indiv, aes(x=delta_cor_obs_median, y=atlas, fill=atlas, colour = atlas))+
#   geom_density_ridges(quantile_lines = F, quantiles = 2, vline_width=0.5, vline_linetype=2, 
#                       alpha=1)+
#   geom_vline(xintercept = 0, lty=1, col="grey80", lwd=0.2)+
#   scale_fill_startrek(guide="none")+
#   scale_color_startrek(guide="none")+
#   coord_cartesian(ylim=c(1.5,4.9))+
#   labs(y="", x="species median Spearman's \U03C1")
#ggsave(paste0("figures/", "cor_obs_ridges_species_median.png"), width=4.5, height = 2.5, bg="white")


ggplot(indiv2, aes(x=delta_cor_obs, y=atlas, fill=species))+
  geom_density_ridges(quantile_lines = F, quantiles = 2, vline_width=0.1, vline_linetype=1, 
                    alpha=0.1, vline_colour = "black", linewidth = 0.1, scale=1)+
  geom_vline(xintercept = 0, lty=1, col="grey80", lwd=0.2)+
  coord_cartesian(ylim=c(1.5,4.5))+
  #scale_fill_startrek(guide="none")+
  labs(y="", x="\U0394 Spearman's \U03C1")+
  theme(legend.position = "none")
ggsave(paste0("figures/", "cor_obs_ridges_species.png"), width=5, height = 4, bg="white", dpi=300)



## jaccard plot ------
jac = readRDS("data/jaccard.rds")
jac$atlas = jac$datasetID
jac$atlas = gsub("26", "Europe", jac$atlas)
jac$atlas = gsub("5", "Czechia", jac$atlas)
jac$atlas = gsub("6", "New York", jac$atlas)
jac$atlas = gsub("17", "New Zealand", jac$atlas)

(jacfig = ggplot(jac, aes(x=jaccard))+
  geom_histogram(fill="grey80")+
  facet_wrap('atlas')+
  scale_fill_startrek(guide="none")+
  scale_color_startrek(guide="none")+
  stat_central_tendency(type="median", aes(col=atlas), lty=2)+
  labs(x='Jaccard between sampling periods, per species', y="count"))
ggsave(paste0("figures/", "jaccard_per_species.png"), width=90, height = 70, units = "mm", bg="white", dpi=300)

(richfig = ggplot(jac, aes(x=friendsTend-friendsT1))+
  geom_histogram(fill="grey80")+ # aes(y = ..density..), 
  facet_grid(time_bin~atlas, scales="free_y")+
#  geom_vline(xintercept = 0, lty=1, col="grey40", lwd=0.3)+
  stat_central_tendency(type="median", aes(col=atlas), lty=2)+
  scale_fill_startrek(guide="none")+
  scale_color_startrek(guide="none")+
  labs(x="Change in number of co-occurring species", y="count")+
  facet_wrap(~atlas, scales="free"))
ggsave(paste0("figures/", "richness_change_per_species.png"), width=70, height = 70, units = "mm", bg="white", dpi=300)



# compare jaccard to mean change per species
indiv = merge(indiv, jac, all.x=TRUE)

ggplot(indiv, aes(x=delta_cor_obs_median, y=jaccard, col=atlas))+
  geom_point(alpha=0.5)+
  facet_wrap('atlas')+
  geom_smooth(col="grey40")+
  labs(x="species median \U0394 Spearman's \U03C1")+
  scale_y_continuous(transform = "exp")+
  scale_color_startrek(guide="none")
ggsave(paste0("figures/", "jaccard_vs_median_cor_obs.png"), width=100, height = 80, units = "mm", bg="white", dpi=150)






### bray curtis -----
bc = readRDS("data/braycurtis.rds")
bc$atlas = gsub("26", "Europe", bc$atlas)
bc$atlas = gsub("5", "Czechia", bc$atlas)
bc$atlas = gsub("6", "New York", bc$atlas)
bc$atlas = gsub("17", "New Zealand", bc$atlas)

(bcfig = ggplot(bc, aes(x=bray_curtis))+
    geom_histogram(fill="grey80")+
    facet_wrap('atlas')+
    scale_fill_startrek(guide="none")+
    scale_color_startrek(guide="none")+
    stat_central_tendency(type="median", aes(col=atlas), lty=2)+
    labs(x='Bray-Curtis between sampling periods, per species', y="count"))
ggsave(paste0("figures/", "bray_curtis_per_species.png"), width=90, height = 70, units = "mm", bg="white", dpi=300)




# add into indiv
indiv = merge(indiv, bc, all.x=TRUE)

ggplot(indiv, aes(x=delta_cor_obs_median, y=1-bray_curtis))+
  geom_point(alpha=0.1)+
  facet_wrap("atlas")








# add avonet -------
avonet = fread("../AVONet/16586228/AVONET1_BirdLife.csv")
names(avonet)
avonet$species = gsub(" ", "_", avonet$Species1)
indiv = merge(indiv, 
              avonet[,.SD, .SDcols = !c('Mass.Refs.Other', 'Reference.species', 'Traits.inferred')],
              by="species", all.x=TRUE)


## function for permutation and plot position
permdiff <- function(x, group, rep = 1000, metric) {
  # x      = data.table
  # group  = column name (string) for grouping variable
  # metric = column name (string) for the metric variable
  
  stopifnot(is.data.table(x))
  stopifnot(group %in% names(x))
  stopifnot("atlas" %in% names(x))
  stopifnot(metric %in% names(x))
  
  # 1) subset to groups with at least 3 members
  tmp <- table(x[[group]])
  valid_groups <- names(tmp[tmp >= 3])
  subx <- x[get(group) %in% valid_groups]
  
  # 2) permutation tests
  group_vec <- c()
  atlas_vec <- c()
  pvals <- c()
  group_median <- c()
  
  for (a in unique(subx$atlas)) {
    sub_atlas <- subx[atlas == a]
    for (g in unique(sub_atlas[[group]])) {
      subset_group <- sub_atlas[get(group) == g]
      if (nrow(subset_group) >= 3) {
        # observed statistic
        obs_median <- median(subset_group[[metric]], na.rm = TRUE)
        # permutation distribution
        perm_stats <- replicate(rep, {
          median(sample(sub_atlas[[metric]], nrow(subset_group)), na.rm = TRUE)
        })
        # two-sided p-value
        perm_pvalue <- mean(abs(perm_stats - median(sub_atlas[[metric]], na.rm=TRUE)) >= 
                              abs(obs_median - median(sub_atlas[[metric]], na.rm=TRUE)))
        
        # collect
        group_vec <- c(group_vec, g)
        atlas_vec <- c(atlas_vec, a)
        pvals <- c(pvals, perm_pvalue)
        group_median <- c(group_median, median(sub_atlas[[metric]]))
      }
    }
  }
  
  # 3) assemble results
  stat_res <- data.table(
    group = group_vec,
    group_median = group_median,
    atlas = atlas_vec,
    p_value = pvals
  )
  
  # safer version for data.table
  stat_res[, significance := fifelse(.SD$p_value < 0.001, "***",
                                     fifelse(.SD$p_value < 0.01, "**",
                                             fifelse(.SD$p_value < 0.05, "*", "")))]
  
  
  # 4) add y-position for plotting
  y_pos <- numeric(nrow(stat_res))
  for (i in seq_len(nrow(stat_res))) {
    ord <- stat_res$group[i]
    atl <- stat_res$atlas[i]
    sub_data <- subx[get(group) == ord & atlas == atl]
    if (nrow(sub_data) > 0) {
      y_pos[i] <- max(sub_data[[metric]], na.rm = TRUE) + 0.02
    } else {
      y_pos[i] <- NA_real_
    }
  }
  stat_res$y <- y_pos
  return(list(data = subx, stats = stat_res))
}


## Taxonomy -----

stat_res = permdiff(indiv, group="Order1", rep=1000, metric="delta_cor_obs_median")
fwrite(stat_res$stats, "output/taxonomy_permutation_stats.csv")

taxplot = ggplot(stat_res$data, aes(y = Order1, x = delta_cor_obs_median, fill = atlas)) +
  geom_vline(data=stat_res$stats, aes(xintercept = group_median), col = "grey") +
  geom_boxplot(varwidth = T, outlier.alpha = 0.2) +
  scale_fill_startrek() +
  labs(y = "", x = "species median \u0394 Spearman's \u03C1") +
  theme(legend.position = "none") +
  facet_wrap(~atlas, ncol = 4, scales="free_x") +
  geom_text(data = stat_res$stats,
            aes(x = y, y = group, label = significance),
            inherit.aes = FALSE,
            hjust = 0,
            size = 4)
ggsave(paste0("figures/", "delta_cor_obs_per_order.png"), width=200, height = 65, unit="mm", bg="white", dpi=300)



## Habitat -----------
stat_res = permdiff(indiv, group="Habitat", rep=1000, metric="delta_cor_obs_median")
fwrite(stat_res$stats, "output/habitat_permutation_stats.csv")

habplot = ggplot(stat_res$data, aes(y = Habitat, x = delta_cor_obs_median, fill = atlas)) +
  geom_vline(data=stat_res$stats, aes(xintercept = group_median), col = "grey") +
  geom_boxplot(varwidth = F, outlier.alpha = 0.2) +
  scale_fill_startrek() +
  labs(y = "", x = "species median \u0394 Spearman's \u03C1") +
  theme(legend.position = "none") +
  facet_wrap(~atlas, ncol = 4, scales="free_x") +
  geom_text(data = stat_res$stats,
            aes(x = y, y = group, label = significance),
            inherit.aes = FALSE,
            hjust = 0,
            size = 4)
ggsave(paste0("figures/", "delta_cor_obs_per_habitat.png"), width=200, height = 40, unit="mm", bg="white", dpi=300)

# ggplot(stat_res$data, aes(y = atlas, x = delta_cor_obs_median, fill = atlas)) +
#   geom_vline(data=stat_res$stats, aes(xintercept = group_median), col = "grey") +
#   geom_boxplot(varwidth = F, outlier.alpha = 0.2) +
#   scale_fill_startrek() +
#   labs(y = "", x = "species median \u0394 Spearman's \u03C1") +
#   theme(legend.position = "none") +
#   facet_wrap(~Habitat, ncol = 4, scales="free_x") #+
#   # geom_text(data = stat_res$stats,
#   #           aes(x = y, y = group, label = significance),
#   #           inherit.aes = FALSE,
#   #           hjust = 0,
#   #           size = 4)
# ggsave(paste0("figures/", "delta_cor_obs_per_habitat_ALT.png"), width=200, height = 40, unit="mm", bg="white", dpi=300)

table(indiv$Habitat, indiv$atlas)


## Primary.Lifestyle ------------------


stat_res = permdiff(indiv, group = "Primary.Lifestyle", rep=10000, metric="delta_cor_obs_median")
fwrite(stat_res$stats, "output/lifestyle_permutation_stats.csv")

ggplot(stat_res$data, aes(y = Primary.Lifestyle, x = delta_cor_obs_median, fill = atlas)) +
  geom_vline(data=stat_res$stats, aes(xintercept = group_median), col = "grey") +
  #geom_boxplot(varwidth = F, outlier.alpha = 0.2) +
  geom_boxplot(varwidth = T, outlier.alpha = 0.2) +
  scale_fill_startrek() +
  labs(y = "", x = "species median \u0394 Spearman's \u03C1") +
  theme(legend.position = "none") +
  facet_grid(~atlas, scales="free") +
  geom_text(data = stat_res$stats, aes(x = y, y = group, label = significance),
            hjust = 0, size = 4)
#ggsave(paste0("figures/", "delta_cor_obs_per_Primary.Lifestyle.png"), width=120, height = 70, unit="mm", bg="white", dpi=300)
ggsave(paste0("figures/", "delta_cor_obs_per_Primary.Lifestyle.png"), width=180, height = 40, unit="mm", bg="white", dpi=300)





## Trophic Niche ------------------
stat_res = permdiff(indiv, group = "Trophic.Niche", rep=10000, metric="delta_cor_obs_median")
fwrite(stat_res$stats, "output/trophicniche_permutation_stats.csv")

trophplot = ggplot(stat_res$data, aes(y = Trophic.Niche, x = delta_cor_obs_median, fill = atlas)) +
  geom_vline(data=stat_res$stats, aes(xintercept = group_median), col = "grey") +
  geom_boxplot(varwidth = F, outlier.alpha = 0.2) +
  scale_fill_startrek() +
  labs(y = "", x = "species median \u0394 Spearman's \u03C1") +
  theme(legend.position = "none") +
  facet_grid(~atlas, scales="free") +
  geom_text(data = stat_res$stats, aes(x = y, y = group, label = significance),
            hjust = 0, size = 4)

ggsave(paste0("figures/", "delta_cor_obs_per_Trophic.Niche.png"), width=200, height = 45, unit="mm", bg="white", dpi=300)


### stats ----
atli = unique(indiv$atlas)
wt = list()
pwt = list()
for(i in 1:length(atli)){
  for(j in 1:length(vars)){}
  tmp = tapply(indiv$delta_cor_obs_median[indiv$atlas==atli[i]], indiv$Trophic.Niche[indiv$atlas==atli[i]], wilcox.test)
  pwt[[i]] = pairwise.wilcox.test(indiv$delta_cor_obs_median[indiv$atlas==atli[i]], indiv$Trophic.Niche[indiv$atlas==atli[i]], p.adjust.method = "fdr")
  tmp2 = rbindlist(tmp)
  tmp2$data.name = names(tmp)
  tmp2$sig = tmp2$p.value<0.05
  tmp2$atlas = atli[i]
  wt[[i]] = tmp2
}
names(pwt) = atli
sig_trophic_niche = list(rbindlist(wt), pwt)
fwrite(sig_trophic_niche[[1]], "output/trophic_niche_stats.csv")


### PLOT
plot_grid(taxplot, habplot, trophplot, ncol=1, labels = c("a)", "b)", "c)"),  
          label_fontface = "plain", label_size = 10, rel_heights = c(65, 40, 45))
ggsave(paste0("figures/", "figS5.png"), width=200, height = 65+40+45, unit="mm", bg="white", dpi=300)


## OTHER -----------------------
ggplot(na.omit(indiv), aes(x=Mass, y=delta_cor_obs_median, col=atlas))+
  geom_point()+scale_x_log10()+geom_smooth()

ggplot(na.omit(indiv), aes(x=Beak.Width, y=delta_cor_obs_median, col=atlas))+
  geom_point()+scale_x_log10()+geom_smooth()

ggplot(na.omit(indiv), aes(x=Beak.Depth, y=delta_cor_obs_median, col=atlas))+
  geom_point()+scale_x_log10()+geom_smooth()

ggplot(na.omit(indiv), aes(x=Tarsus.Length, y=delta_cor_obs_median, col=atlas))+
  geom_point()+scale_x_log10()+geom_smooth()

ggplot(na.omit(indiv), aes(x=Wing.Length, y=delta_cor_obs_median, col=atlas))+
  geom_point()+scale_x_log10()+geom_smooth()

ggplot(na.omit(indiv), aes(x=Centroid.Latitude, y=delta_cor_obs_median, col=atlas))+
  geom_point()+scale_x_log10()+geom_smooth()

ggplot(na.omit(indiv), aes(x=Range.Size, y=delta_cor_obs_median, col=atlas))+
  geom_point()+scale_x_log10()+geom_smooth()








### no median, all values ---------------------------

# this uses the data.table before getting medians per species. 
indiv2 = merge(indiv2, 
              avonet[,.SD, .SDcols = !c('Mass.Refs.Other', 'Reference.species', 'Traits.inferred')],
              by="species", all.x=TRUE)

ggplot(na.omit(indiv2), aes(x=delta_cor_obs, y=atlas, fill=atlas, col=atlas))+
  geom_density_ridges(quantile_lines = F, quantiles = 2, vline_width=0.1, vline_linetype=1, 
                      alpha=1, vline_colour = "black", linewidth = 0.1)+
  geom_vline(xintercept = 0, lty=1, col="grey80", lwd=0.2)+
  coord_cartesian(ylim=c(1.5,5.2))+
  scale_fill_startrek(guide="none")+
  scale_color_startrek(guide="none")+
  labs(y="", x="\U0394 Spearman's \U03C1")+
  theme(legend.position = "bottom")
# looks the same as chan - as it should

ggplot(na.omit(indiv2), aes(x=delta_cor_obs, y=Primary.Lifestyle, fill=atlas))+
  geom_vline(xintercept = 0, col = "grey") +
  geom_boxplot(varwidth = F, outlier.alpha = 0.2) +
  scale_fill_startrek() +
  labs(y = "", x = "species median \u0394 Spearman's \u03C1") +
  theme(legend.position = "none") +
  facet_wrap(~atlas, ncol = 4)

atli = unique(indiv2$atlas)
wt = list()
pwt = list()
for(i in 1:length(atli)){
  for(j in 1:length(vars)){}
  tmp = tapply(indiv2$delta_cor_obs[indiv2$atlas==atli[i]], indiv2$Primary.Lifestyle[indiv2$atlas==atli[i]], wilcox.test)
  pwt[[i]] = pairwise.wilcox.test(indiv2$delta_cor_obs[indiv2$atlas==atli[i]], indiv2$Primary.Lifestyle[indiv2$atlas==atli[i]], p.adjust.method = "fdr")
  tmp2 = rbindlist(tmp)
  tmp2$data.name = names(tmp)
  tmp2$sig = tmp2$p.value<0.05
  tmp2$atlas = atli[i]
  wt[[i]] = tmp2
}
names(pwt) = atli
sig_lifestyle = list(rbindlist(wt), pwt)
sig_lifestyle


### The 2nd species ------
indiv2$sp2 <- mapply(function(pair, sp) gsub(sp, "", pair, fixed = TRUE),
                     indiv2$species_pair, indiv2$species)
# do not do a simple gsub regex grab, order might be off (it should not but you neer know)
indiv2$sp2 = gsub("\\|", "", indiv2$sp2)
indiv2 = merge(indiv2, 
               avonet[,.SD, .SDcols = !c('Mass.Refs.Other', 'Reference.species', 'Traits.inferred')],
               by.x="sp2", by.y="species", all.x=TRUE)

ggplot(na.omit(indiv2), aes(x=delta_cor_obs, y=Primary.Lifestyle.x, fill=Primary.Lifestyle.y))+
  geom_vline(xintercept = 0, lty=1, col="grey80", lwd=0.5)+
  geom_boxplot(varwidth = F, outlier.alpha = 0.05)+
  scale_fill_startrek(name="2^nd^ species")+
  scale_color_startrek(guide="none")+
  labs(y="1^st^ species", x="\U0394 Spearman's \U03C1")+
  facet_wrap("atlas")+
  theme(legend.position = "bottom")
ggsave(paste0("figures/", "delta_cor_obs_per_ecology.png"), width=200, height = 150, unit="mm", bg="white", dpi=150)


### stats -------------
atli = unique(indiv2$atlas)
wt = list()
for(i in 1:length(atli)){
  tmp = indiv2[indiv2$atlas==atli[i],]
  sp1_lifestyle = sort(unique(tmp$Primary.Lifestyle.x))
  
  for(j in 1:length(sp1_lifestyle)){
    tmp_wt = tapply(tmp$delta_cor_obs[tmp$Primary.Lifestyle.x==sp1_lifestyle[j]], 
                  tmp$Primary.Lifestyle.y[tmp$Primary.Lifestyle.x==sp1_lifestyle[j]], wilcox.test)
    tmp2 = rbindlist(tmp_wt)
    tmp2$data.name = names(tmp_wt)
    tmp2$sp1_primary_lifestyle = sp1_lifestyle[j]
    if(j==1){wt_res = tmp2}else{wt_res = rbind(wt_res, tmp2)}
  }
  wt_res$atlas = atli[i]
  wt_res$sig = wt_res$p.value<0.05
  wt[[i]] = wt_res
}
sig_ecology = rbindlist(wt)
sig_lifestyle


tmp1 = pairwise.wilcox.test(indiv2$delta_cor_obs[indiv2$atlas=="Czechia" & indiv2$Primary.Lifestyle.x=="Aquatic"], 
                     indiv2$Primary.Lifestyle.y[indiv2$atlas=="Czechia" & indiv2$Primary.Lifestyle.x=="Aquatic"], p.adjust.method = "fdr")
tmp2 = pairwise.wilcox.test(indiv2$delta_cor_obs[indiv2$atlas=="Czechia" & indiv2$Primary.Lifestyle.x=="Terrestrial"], 
                     indiv2$Primary.Lifestyle.y[indiv2$atlas=="Czechia" & indiv2$Primary.Lifestyle.x=="Terrestrial"], p.adjust.method = "fdr")
tmp2 = pairwise.wilcox.test(indiv2$delta_cor_obs[indiv2$atlas=="Czechia" & indiv2$Primary.Lifestyle.x=="Terrestrial"], 
                            indiv2$Primary.Lifestyle.y[indiv2$atlas=="Czechia" & indiv2$Primary.Lifestyle.x=="Terrestrial"], p.adjust.method = "fdr")

# ggplot(na.omit(indiv2[Primary.Lifestyle.x=="Aquatic"]), aes(x=delta_cor_obs, y=Primary.Lifestyle.x, fill=Primary.Lifestyle.y))+
#   geom_vline(xintercept = 0, lty=1, col="grey80", lwd=0.5)+
#   geom_boxplot(varwidth = T, outlier.alpha = 0.05)+
#   geom_vline(xintercept = 0, lty=1, col="grey60", lwd=0.5)+
#   scale_fill_viridis_d("")+
#   coord_cartesian(expand = FALSE)+
#   labs(y="", x="\U0394 Spearman's \U03C1")+
#   facet_wrap("atlas")+
#   theme(legend.position = "bottom")

# ggplot(na.omit(indiv2), aes(x=delta_cor_obs, y=Primary.Lifestyle.x, fill=Primary.Lifestyle.y, col=Primary.Lifestyle.y))+
#   geom_density_ridges(quantile_lines = F, quantiles = 2, vline_width=0.1, vline_linetype=1, 
#                       alpha=0.3, vline_colour = "black", linewidth = 0.1)+
#   geom_vline(xintercept = 0, lty=1, col="grey80", lwd=0.2)+
#   #coord_cartesian(ylim=c(1.5,5.3))+
# #  scale_fill_manual(values=startrek, guide="none")+
# #  scale_color_manual(values=startrek, guide="none")+
#   labs(y="", x="\U0394 Spearman's \U03C1")+
#   facet_wrap("atlas")+
#   theme(legend.position = "bottom")

# ggplot(na.omit(indiv2), aes(x=Mass.x, y=Primary.Lifestyle.x, fill=atlas))+
#   geom_boxplot(varwidth = T)+
#   geom_vline(xintercept = 0, lty=1, col="grey80", lwd=0.2)+
#   coord_cartesian(ylim=c(1.5,5.3))+
#   scale_fill_startrek(guide="none")+
#   scale_color_startrek(guide="none")+
#   labs(y="", x="\U0394 Spearman's \U03C1")+
#   theme(legend.position = "bottom")
# 


















# PHYLOGENY -------------------------------------------------------------------------
# for czechia
# indivCZ = indiv[atlas=="Czechia",]
# phyCZ = keep.tip(phy, tip=unique(indivCZ$species)[unique(indivCZ$species) %in% phy$tip.label])
# phyCZ$median_correlation_change = indivCZ$delta_cor_obs_median[match(phyCZ$tip.label, indivCZ$species)]
# phyCZ$Order = indivCZ$Order1[match(phyCZ$tip.label, indivCZ$species)]
# phyorder = as.numeric(as.factor(phyCZ$Order))
# distinct_colors <- c(
#   "#E41A1C",  # Red
#   "#377EB8",  # Blue
#   "#4DAF4A",  # Green
#   "#984EA3",  # Purple
#   "#FF7F00",  # Orange
#   "#FFFF33",  # Yellow
#   "#A65628",  # Brown
#   "#F781BF",  # Pink
#   "#999999",  # Gray
#   "#66C2A5",  # Teal
#   "#FC8D62",  # Salmon
#   "#8DA0CB",  # Light blue
#   "#E78AC3",  # Light pink
#   "#A6D854",  # Lime green
#   "#FFD92F",  # Bright yellow
#   "#E5C494"   # Tan
# )
# 
# 
# plot.phylo(phyCZ, type="fan", 
#            tip.color = distinct_colors[phyorder], cex=0.5)
# 
# red_to_blue <- colorRampPalette(c("blue", "red"))(6)
# medchange = cut(phyCZ$median_correlation_change, 6)
# plot.phylo(phyCZ, type="fan", 
#            tip.color = red_to_blue[as.numeric(medchange)], cex=0.5)
# 
# plot.phylo(phyCZ, type="fan", 
#            tip.color = c("red", "blue")[as.numeric(phyCZ$median_correlation_change<=0)+1], cex=0.5)








## Range --------------------------------------------


# merge change in occupancy into indiv data
tmp = unique(dat[, .(scientificName, change_occupancy, change_pres.abs, datasetID)])
setnames(tmp, "scientificName", "species")
tmp$atlas = tmp$datasetID
tmp$atlas = gsub(5, "Czechia", tmp$atlas)
tmp$atlas = gsub("26", "Europe", tmp$atlas)
tmp$atlas = gsub("6", "New York", tmp$atlas)
tmp$atlas = gsub("17", "New Zealand", tmp$atlas)

indiv = merge(indiv, tmp[,.(species, atlas, change_occupancy, change_pres.abs)], by=c("atlas", "species"), 
              all.x=TRUE)

# log ratio for delta occupancy
library(compositions)
clr_data <- clr(indiv$change_occupancy) 

(occ1 = ggplot(indiv, aes(x=as.numeric(clr(change_occupancy)), y=delta_cor_obs_median, col=atlas, fill=atlas))+
  geom_vline(xintercept = 0, col='grey') + geom_hline(yintercept = 0, col='grey')+
  geom_point(alpha=0.5)+
  scale_color_startrek()+
  scale_fill_startrek()+
  geom_smooth(alpha=0.1, lwd=0.5)+
#  stat_cor(method = "spearman", size=3, p.accuracy = 0.001, cor.coef.name = "rho")+
  labs(x="\U0394 occupancy (centered log-ratio)", y="median \U0394 Spearman's \U03C1")+
  theme(legend.position = "bottom"))
(occ2 = ggplot(indiv, aes(x=as.numeric(clr(change_occupancy)), y=delta_cor_obs_median, col=atlas))+
    geom_vline(xintercept = 0, col='grey') + geom_hline(yintercept = 0, col='grey')+
    geom_point(alpha=0.7)+
    scale_color_startrek(guide="none")+
    geom_smooth(alpha=0.5, lwd=0.5)+
    stat_cor(method = "spearman", size=3, p.accuracy = 0.001, cor.coef.name = "rho",
             label.sep = "\n", col="black")+
    facet_wrap(~atlas, ncol=2)+
    labs(x="\U0394 occupancy (centered log-ratio)", y="median \U0394 Spearman's \U03C1")+
    theme(legend.position = "bottom"))
ggsave(paste0("figures/", "delta_occupancy_delta_cor_lifestyle.png"), units="mm", width=120, height = 100, bg="white", dpi = 300)


cor.test(indiv$change_occupancy, indiv$delta_cor_obs_median, method="s")
tapply(indiv, indiv$atlas, function(x){cor.test(x$delta_cor_obs_median, as.numeric(clr(x$change_occupancy)), method='s')})





# merge in range size (n occupied cells)
dat[, grid_n:=sum(pres.abs), by=.(datasetID, scalingID, endYear, scientificName)]
tmp = unique(dat[, .(grid_n, datasetID, scalingID, endYear, scientificName, mean.psi_sum)])
tmp$atlas = tmp$datasetID
tmp$atlas = gsub(5, "Czechia", tmp$atlas)
tmp$atlas = gsub("26", "Europe", tmp$atlas)
tmp$atlas = gsub("6", "New York", tmp$atlas)
tmp$atlas = gsub("17", "New Zealand", tmp$atlas)
setnames(tmp, "scientificName", "species")
indiv = merge(indiv, tmp[endYear=="1", .(grid_n, atlas, species)], all.x=TRUE)

ggplot(indiv, aes(y=Primary.Lifestyle, x=grid_n, fill=atlas))+
  geom_boxplot(varwidth = F)+
  scale_fill_startrek(guide="none")+
  labs(y="", x="N grids")+
  facet_wrap(~atlas, ncol=4)
ggsave(paste0("figures/", "Ngrids_lifestyle.png"), units="mm", width=180, height = 40, bg="white", dpi = 300)

tapply(indiv, indiv$atlas, function(x){pairwise.wilcox.test(x$grid_n, x$Primary.Lifestyle, p.adjust.method = "fdr")})
pairwise.wilcox.test(indiv$grid_n, indiv$Primary.Lifestyle, p.adjust.method = "fdr")







### occupancy change per order --------------
stat_res = permdiff(indiv, group="Order1", rep=1000, metric="change_occupancy")
fwrite(stat_res$stats, "output/occupancy_taxonomy_permutation_stats.csv")

#ggplot(stat_res$data, aes(y = Order1, x = as.numeric(clr(change_occupancy)), fill = atlas)) +
ggplot(stat_res$data, aes(y = Order1, x = change_occupancy, fill = atlas)) +
  geom_vline(data=stat_res$stats, aes(xintercept = group_median), col = "grey") +
  geom_boxplot(varwidth = F, outlier.alpha = 0.2) +
  scale_fill_startrek() +
  labs(y = "", x = "\U0394 occupancy relative to 1^st^ sampling period") +
  theme(legend.position = "none") +
  facet_wrap(~atlas, ncol = 4, scales="free_x") +
  geom_text(data = stat_res$stats,
            aes(x = y, y = group, label = significance),
            inherit.aes = FALSE,
            hjust = 0,
            size = 4)

ggsave(paste0("figures/", "occupancy_change_per_order.png"), width=200, height = 65, unit="mm", bg="white", dpi=150)






### occupancy change per habitat --------------
stat_res = permdiff(indiv, group="Habitat", rep=1000, metric="change_occupancy")
fwrite(stat_res$stats, "output/occupancy_habitat_permutation_stats.csv")

ggplot(stat_res$data, aes(y = Habitat, x = change_occupancy, fill = atlas)) +
  geom_vline(data=stat_res$stats, aes(xintercept = group_median), col = "grey") +
  geom_boxplot(varwidth = F, outlier.alpha = 0.2) +
  scale_fill_startrek() +
  labs(y = "", x = "\U0394 occupancy relative to 1^st^ sampling period (log ratio)") +
  theme(legend.position = "none") +
  facet_wrap(~atlas, ncol = 4, scales="free_x") +
  geom_text(data = stat_res$stats,
            aes(x = y, y = group, label = significance),
            inherit.aes = FALSE,
            hjust = 0,
            size = 4)

ggsave(paste0("figures/", "occupancy_change_per_habitat.png"), width=200, height = 40, unit="mm", bg="white", dpi=150)






### occupancy change per primary lifestyle --------------
stat_res = permdiff(indiv, group="Primary.Lifestyle", rep=1000, metric="change_occupancy")
fwrite(stat_res$stats, "output/occupancy_habitat_permutation_stats.csv")

ggplot(stat_res$data, aes(y = Primary.Lifestyle, x = change_occupancy, fill = atlas)) +
  geom_vline(data=stat_res$stats, aes(xintercept = group_median), col = "grey") +
  geom_boxplot(varwidth = F, outlier.alpha = 0.2) +
  scale_fill_startrek() +
  scale_x_continuous(trans="log1p")+
  labs(y = "", x = "\U0394 occupancy relative to 1^st^ sampling period") +
  theme(legend.position = "none") +
  facet_wrap(~atlas, ncol = 4, scales="free_x") +
  geom_text(data = stat_res$stats,
            aes(x = y, y = group, label = significance),
            inherit.aes = FALSE,
            hjust = 0,
            size = 4)

ggsave(paste0("figures/", "occupancy_change_per_lifestyle.png"), width=180, height = 40, unit="mm", bg="white", dpi=300)














# Balanced ranges -----
# chan[ , balanced_range_ratio:=range_ratio_pres.abs>0.1]
# chan[ , upper_range:=max(c(pres.abs_sum.x, pres.abs_sum.y)), by=atlas]
# chan[, average_ranges:= pres.abs_sum.x>50 & pres.abs_sum.y>50 & # both species above 50
#        pres.abs_sum.x<upper_range*0.9 & pres.abs_sum.y<upper_range*0.9, ]
chan[, small_ranges:= pres.abs_sum.x<50 | pres.abs_sum.y<50, ]

ggplot(chan, aes(x=range_ratio_pres.abs))+geom_histogram()+facet_grid(~atlas)
tapply(chan$range_ratio_pres.abs>0.15, chan$atlas, table) # find a good balance
tapply(chan$range_ratio_pres.abs>0.15, chan$atlas, table) # find a good balance


sub = chan[balanced_range_ratio>=0.15 & small_ranges==FALSE,]
tapply(sub$atlas, sub$sp1, length)

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
    # if(i==1){
    #   indiv_int2 = data.table(delta_cor_obs = tmp$delta_cor_obs[id],
    #                           species = indiv_int$species[i],
    #                           atlas = at[a],
    #                           species_pair = tmp$species_pair[id])
    # }else{tmp2 = data.table(delta_cor_obs = tmp$delta_cor_obs[id],
    #                         species = indiv_int$species[i],
    #                         atlas = at[a],
    #                         species_pair = tmp$species_pair[id])
    # indiv_int2 = rbind(indiv_int2, tmp2)}
    
  }
  indiv_br[[a]] = indiv_int
  #indiv2[[a]] = indiv_int2
  cat(a, i,  "\r")
}
indiv_br = rbindlist(indiv_br)


# add avonet
indiv_br = merge(indiv_br, 
              avonet[,.SD, .SDcols = !c('Mass.Refs.Other', 'Reference.species', 'Traits.inferred')],
              by="species", all.x=TRUE)

# indiv_nosmall = na.omit(indiv[grid_n>50])
stat_res = permdiff(indiv_br, group="Primary.Lifestyle", rep=1000, metric="delta_cor_obs_median")
fwrite(stat_res$stats, "output/primarylifestyle_permutation_stats_nosmall.csv")

ggplot(stat_res$data, aes(y = Primary.Lifestyle, x = delta_cor_obs_median, fill = atlas)) +
  geom_vline(data=stat_res$stats, aes(xintercept = group_median), col = "grey") +
  geom_boxplot(varwidth = T, outlier.alpha = 0.2) +
  scale_fill_startrek() +
  labs(y = "", x = "species median \u0394 Spearman's \u03C1") +
  theme(legend.position = "none") +
  facet_wrap(~atlas, ncol = 4, scales="free_x") +
  geom_text(data = stat_res$stats,
            aes(x = y, y = group, label = significance),
            inherit.aes = FALSE,
            hjust = 0,
            size = 4)

ggsave(paste0("figures/", "Primary.Lifestyle_larger_ranges.png"), width=200, height = 40, unit="mm", bg="white", dpi=300)

#ggsave(paste0("figures/", "balanced_pairs.png"), width=200, height = 60, unit="mm", bg="white", dpi=150)











## species in multiple atlas --------------


tmp <- indiv[species %in% names(which(table(indiv$species)>1)),]

ggplot(tmp, aes(y=species, x=delta_cor_obs_median, col=atlas))+
  geom_point()+
  theme(panel.grid.major.y = element_line())












# example maps -----------------
gc()
dataset_id
grid = lapply(paste0("data/all_scales_atlas_", dataset_id, ".gpkg"), vect)
grid = vect(grid)

# czechoutline = grid[grid$scalingID==64 & ,1]
# czechoutline = st_as_sf(czechoutline)

# species pairs with positive delta cor
res = res[over==TRUE,]
res$datasetID = res$atlas
res$datasetID = gsub("Czechia", "5", res$datasetID)
res$datasetID = gsub("New York", "6", res$datasetID)


# increase / decreases?
ord = TRUE
head(res[order(res$delta_cor_obs, decreasing = ord), .(species_pair, Order1.x, Order1.y, overlap, atlas, delta_cor_obs)], 20)


for(i in 1:100){
  sp = res$species_pair[order(res$delta_cor_obs, decreasing = ord)][i]
  dID = res$datasetID[order(res$delta_cor_obs, decreasing = ord)][i]Ftende
  change = res$delta_cor_obs[order(res$delta_cor_obs, decreasing = ord)][i]
  sp = unlist(strsplit(sp, "|", fixed = TRUE))
  tmp = dat[dat$scientificName %in% sp & dat$datasetID %in% dID,]
  tmp = tmp[,.(pres.abs, mean.psi, scientificName, endYear, siteID, scalingID, datasetID)]
  # grid$datasetID = as.character(grid$datasetID)
  # tmp$scalingID = as.integer(tmp$scalingID)
  grid_sub = grid[grid$datasetID==dID,]
  p1 = merge(grid_sub, tmp, by=c("datasetID", "scalingID", "siteID"), all.x=TRUE)
  p1 = st_as_sf(p1)
  p1$years = p1$endYear
  if(dID == "5"){
    p1$years = gsub("1", "1985-1989", p1$years) 
    p1$years = gsub("3", "2014-2017", p1$years)
    p1 = p1[p1$endYear %in% c("1", "3"),]
  }
  if(dID == "6"){ # 1980-1985, 2000-2005
    p1$years = gsub("1", "1980-1985", p1$years) 
    p1$years = gsub("2", "2000-2005", p1$years)
  }
  
  # separated
  ggplot(na.omit(p1[p1$pres.abs==1,]), 
         aes(fill=factor(pres.abs)))+
    #geom_sf(data=st_as_sf(grid_sub[grid_sub$datasetID==dID & grid_sub$scalingID==max(grid_sub$scalingID),])
    #        , fill=NA, col="grey50")+
    geom_sf(alpha=1)+
    scale_fill_startrek(name="presence")+
    facet_grid(scientificName~years)+
    theme_map()+
    theme(legend.position = "none")+
    ggtitle(paste0("delta cor = ",  round(change, 2)))+
    theme(title = element_text(size=8, face = "plain"))
  ggsave(paste0("figures/maps/", paste(sp, collapse = "-"), ".png"), width=7, height = 4.5, bg="white")
  
  # occupancy
  ggplot(na.omit(p1), 
         aes(fill=mean.psi))+
    #geom_sf(data=st_as_sf(grid_sub[grid_sub$datasetID==dID & grid_sub$scalingID==max(grid_sub$scalingID),])
    #        , fill=NA, col="grey50")+
    geom_sf(alpha=1)+
    scale_fill_viridis_c(option = "F")+
    facet_grid(scientificName~years)+
    theme_map()+
    theme(legend.position = "bottom")+
    ggtitle(paste0("delta cor = ",  round(change, 2)))+
    theme(title = element_text(size=8, face = "plain"), 
          legend.text = element_text(size=8))
  ggsave(paste0("figures/maps/", paste(sp, collapse = "-"), "_conti.png"), width=7, height = 4.5, bg="white")
  cat(i, "\r")
}














##############################################################################




ggplot(na.omit(indiv), aes(x=grid_n, y=delta_cor_obs_median, col=Primary.Lifestyle))+
  geom_point(size=2, alpha=0.5)+
  geom_hline(yintercept = 0, col="grey20", lty=2)+
  stat_smooth(se = FALSE)+
  scale_color_startrek()+
  labs(x="presence grids T1", y="median \U0394 correlation per species")+
  facet_wrap("atlas", scales="free_x")+
  theme(legend.position = "bottom")

ggplot(na.omit(indiv), aes(y=Primary.Lifestyle, x=grid_n, fill=Primary.Lifestyle))+
  geom_boxplot(varwidth = T, outlier.alpha = 0.2)+
  scale_fill_startrek()+
  labs(y="", x="presence grids T1")+
  facet_wrap("atlas")


# range size vs range size ratio
indiv_range = list()
at = sort(unique(res$atlas))
for(a in 1:length(at)){
  tmp = res[atlas==at[a],]
  indiv_int = data.table(species = sort(unique(c(tmp$sp1, tmp$sp2))), 
                         range_ratio_median = NA, 
                         atlas = at[a])
  for(i in 1:nrow(indiv_int)){
    id = grep(indiv_int$species[i], tmp$species_pair)
    indiv_int$range_ratio_median[i] = median(tmp$range_ratio_pres.abs[id], na.rm=TRUE)
  }
  indiv_range[[a]] = indiv_int
  cat(a, i,  "\r")
}
indiv_range = rbindlist(indiv_range)
indiv = merge(indiv, indiv_range, all.x=TRUE)

ggplot(indiv, aes(y=delta_cor_obs_median, x=range_ratio_median, col=grid_n))+
  geom_point(alpha=0.4)+
  facet_wrap("atlas", scales="free_x")+
  theme(legend.position = "bottom")

ggplot(indiv, aes(x=grid_n, y=range_ratio_median, col=atlas))+
  geom_point()

tapply(indiv$range_ratio_median, indiv$atlas, median, na.rm=TRUE)


ggplot(res[pres.abs_sum.x>100 & pres.abs_sum.y>100 &
             pres.abs_sum.x<500 & pres.abs_sum.y<500,], aes(x=delta_cor_obs))+
  geom_histogram()+
  facet_wrap("atlas", scales="free_x")





inp = res[pres.abs_sum.x>50 & pres.abs_sum.y>50 &
            pres.abs_sum.x<500 & pres.abs_sum.y<500,]
ord = FALSE
# maps 
for(i in 1:100){
  sp = inp$species_pair[order(inp$delta_cor_obs, decreasing = ord)][i]
  dID = inp$datasetID[order(inp$delta_cor_obs, decreasing = ord)][i]
  change = inp$delta_cor_obs[order(inp$delta_cor_obs, decreasing = ord)][i]
  sp = unlist(strsplit(sp, "|", fixed = TRUE))
  tmp = dat[dat$scientificName %in% sp & dat$datasetID %in% dID,]
  tmp = tmp[,.(pres.abs, mean.psi, scientificName, endYear, siteID, scalingID, datasetID)]
  # grid$datasetID = as.character(grid$datasetID)
  # tmp$scalingID = as.integer(tmp$scalingID)
  grid_sub = grid[grid$datasetID==dID,]
  p1 = merge(grid_sub, tmp, by=c("datasetID", "scalingID", "siteID"), all.x=TRUE)
  p1 = st_as_sf(p1)
  p1$years = p1$endYear
  if(dID == "5"){
    p1$years = gsub("1", "1985-1989", p1$years) 
    p1$years = gsub("3", "2014-2017", p1$years)
    p1 = p1[p1$endYear %in% c("1", "3"),]
  }
  if(dID == "6"){ # 1980-1985, 2000-2005
    p1$years = gsub("1", "1980-1985", p1$years) 
    p1$years = gsub("2", "2000-2005", p1$years)
  }
  
  # separated
  ggplot(na.omit(p1[p1$pres.abs==1,]), 
         aes(fill=factor(pres.abs)))+
    #geom_sf(data=st_as_sf(grid_sub[grid_sub$datasetID==dID & grid_sub$scalingID==max(grid_sub$scalingID),])
    #        , fill=NA, col="grey50")+
    geom_sf(alpha=1)+
    scale_fill_startrek(name="presence")+
    facet_grid(scientificName~years)+
    theme_map()+
    theme(legend.position = "none")+
    ggtitle(paste0("delta cor = ",  round(change, 2)))+
    theme(title = element_text(size=8, face = "plain"))
  ggsave(paste0("figures/maps/", paste(sp, collapse = "-"), "_balanced_ranges.png"), width=7, height = 4.5, bg="white")
  
  # occupancy
  ggplot(na.omit(p1), 
         aes(fill=mean.psi))+
    #geom_sf(data=st_as_sf(grid_sub[grid_sub$datasetID==dID & grid_sub$scalingID==max(grid_sub$scalingID),])
    #        , fill=NA, col="grey50")+
    geom_sf(alpha=1)+
    scale_fill_viridis_c(option = "F")+
    facet_grid(scientificName~years)+
    theme_map()+
    theme(legend.position = "bottom")+
    ggtitle(paste0("delta cor = ",  round(change, 2)))+
    theme(title = element_text(size=8, face = "plain"), 
          legend.text = element_text(size=8))
  ggsave(paste0("figures/maps/", paste(sp, collapse = "-"), "_conti_balanced_ranges.png"), width=7, height = 4.5, bg="white")
  cat(i, "\r")
}




























### PLOT THE AQUATICS -----------------------------------------------------



sel = indiv$species[indiv$Primary.Lifestyle %in% c("Aquatic")]
inp = res[grep(paste0(sel, collapse = "|"), res$species_pair), ]

ord = TRUE
# maps 
for(i in 1:100){
  sp = inp$species_pair[order(inp$delta_cor_obs, decreasing = ord)][i]
  dID = inp$datasetID[order(inp$delta_cor_obs, decreasing = ord)][i]
  change = inp$delta_cor_obs[order(inp$delta_cor_obs, decreasing = ord)][i]
  sp = unlist(strsplit(sp, "|", fixed = TRUE))
  tmp = dat[dat$scientificName %in% sp & dat$datasetID %in% dID,]
  tmp = tmp[,.(pres.abs, mean.psi, scientificName, endYear, siteID, scalingID, datasetID)]
  # grid$datasetID = as.character(grid$datasetID)
  # tmp$scalingID = as.integer(tmp$scalingID)
  grid_sub = grid[grid$datasetID==dID,]
  p1 = merge(grid_sub, tmp, by=c("datasetID", "scalingID", "siteID"), all.x=TRUE)
  p1 = st_as_sf(p1)
  p1$years = p1$endYear
  if(dID == "5"){
    p1$years = gsub("1", "1985-1989", p1$years) 
    p1$years = gsub("3", "2014-2017", p1$years)
    p1 = p1[p1$endYear %in% c("1", "3"),]
  }
  if(dID == "6"){ # 1980-1985, 2000-2005
    p1$years = gsub("1", "1980-1985", p1$years) 
    p1$years = gsub("2", "2000-2005", p1$years)
  }
  
  # separated
  ggplot(na.omit(p1[p1$pres.abs==1,]), 
         aes(fill=factor(pres.abs)))+
    #geom_sf(data=st_as_sf(grid_sub[grid_sub$datasetID==dID & grid_sub$scalingID==max(grid_sub$scalingID),])
    #        , fill=NA, col="grey50")+
    geom_sf(alpha=1)+
    scale_fill_startrek(name="presence")+
    facet_grid(scientificName~years)+
    theme_map()+
    theme(legend.position = "none")+
    ggtitle(paste0("delta cor = ",  round(change, 2)))+
    theme(title = element_text(size=8, face = "plain"))
  ggsave(paste0("figures/maps/AQUATICS/", paste(sp, collapse = "-"), ".png"), width=7, height = 4.5, bg="white")
  
  # occupancy
  ggplot(na.omit(p1), 
         aes(fill=mean.psi))+
    #geom_sf(data=st_as_sf(grid_sub[grid_sub$datasetID==dID & grid_sub$scalingID==max(grid_sub$scalingID),])
    #        , fill=NA, col="grey50")+
    geom_sf(alpha=1)+
    scale_fill_viridis_c(option = "F")+
    facet_grid(scientificName~years)+
    theme_map()+
    theme(legend.position = "bottom")+
    ggtitle(paste0("delta cor = ",  round(change, 2)))+
    theme(title = element_text(size=8, face = "plain"), 
          legend.text = element_text(size=8))
  ggsave(paste0("figures/maps/AQUATICS/", paste(sp, collapse = "-"), "_conti.png"), width=7, height = 4.5, bg="white")
  cat(i, "\r")
}



# what order is the second species

aq = chan
tmp = unique(indiv[, .(species, Primary.Lifestyle)])
aq = merge(aq, tmp, by.x='sp1', by.y="species", all.x=TRUE)
aq = merge(aq, tmp, by.x='sp2', by.y="species", all.x=TRUE)


ggplot(aq[Primary.Lifestyle.y=="Aquatic",], aes(x=delta_cor_obs, y=atlas, fill=Primary.Lifestyle.x))+
  geom_boxplot(outlier.alpha = 0.05, varwidth = T)+
  geom_vline(xintercept = 0, lty=1, col="grey60", lwd=0.5)+
  labs(y="", x="\U0394 Spearman's \U03C1")+
  coord_cartesian(ylim = c(1.3,1.7))
#ggsave(paste0("figures/", "delta_cor_obs_aquatics_orders.png"), width=4.5, height = 2.5, bg="white")



other_order = gsub(paste0(sel, collapse = "|"), "", inp$species_pair)
tmp = as.data.table(sort(table(gsub("\\|", "", other_order))))
tmp = merge(tmp, avonet[,.(species, Primary.Lifestyle)], by.x="V1", by.y="species", all.x=TRUE)
tmp$V1[1] = "aquatic"
tmp$Primary.Lifestyle[1] = "aquatic"
tapply(tmp$N, tmp$Primary.Lifestyle, sum)/sum(tmp$N)


# how many per primary lifestyle in the atlases
tapply(indiv$Primary.Lifestyle, indiv$atlas, function(x){table(x)/length(x)})


# range size in general across lifestyles
ggplot(na.omit(indiv), aes(y=Primary.Lifestyle, x=Range.Size))+
  geom_boxplot()
ggplot(na.omit(indiv), aes(y=Primary.Lifestyle, x=grid_n, fill=atlas))+
  geom_boxplot(varwidth = T)+
  scale_fill_startrek()
# tapply(indiv$grid_n[indiv$atlas=="Czechia"], indiv$Primary.Lifestyle[indiv$atlas=="Czechia"], wilcox.test)
# tapply(indiv$grid_n[indiv$atlas=="New York"], indiv$Primary.Lifestyle[indiv$atlas=="New York"], wilcox.test)
pairwise.wilcox.test(indiv$grid_n[indiv$atlas=="Czechia"], indiv$Primary.Lifestyle[indiv$atlas=="Czechia"], p.adjust.method = "fdr")
pairwise.wilcox.test(indiv$grid_n[indiv$atlas=="New York"], indiv$Primary.Lifestyle[indiv$atlas=="New York"], p.adjust.method = "fdr")


# percent of the positive changes
table(indiv$Primary.Lifestyle)/length(indiv$Primary.Lifestyle)

length(indiv$Primary.Lifestyle[indiv$delta_cor_obs_median>0.1])
length(indiv$Primary.Lifestyle[indiv$delta_cor_obs_median<=(-0.1)])

table(indiv$Primary.Lifestyle[indiv$delta_cor_obs_median>0.1])/
  length(indiv$Primary.Lifestyle[indiv$delta_cor_obs_median>0.1])
table(indiv$Primary.Lifestyle[indiv$delta_cor_obs_median<=0.1])/
  length(indiv$Primary.Lifestyle[indiv$delta_cor_obs_median<=0.1])

quantile(indiv$delta_cor_obs_median, probs=c(0.05, 0.5, 0.85, 0.90))


# what are the average c scores per Primary /lifestyle?






















# no overlap one:


### plot Accipiter_gentilis|Accipiter_nisus
dataset_id = 6
dat = fread(paste0("data/all_scales_atlas_", dataset_id, ".csv"))
grid = vect(paste0("data/all_scales_atlas_", dataset_id, ".gpkg"))
dat$scientificName = gsub(pattern = " ", "_", dat$scientificName)

sp1 = "Accipiter_cooperii"
sp2 = "Accipiter_gentilis"
tmp = dat[scientificName %in% c(sp1,sp2),]
tmp = unique(tmp[,.(scientificName, endYear, siteID, datasetID)])
p1 = merge(grid, tmp[,.(scientificName, endYear, siteID, datasetID)], all.x=TRUE)

ggplot(na.omit(st_as_sf(p1[p1$scalingID==1,])), aes(fill=scientificName))+
  geom_sf(alpha=0.4, col=NA)+
  facet_wrap(~endYear)+
  theme(legend.position = "bottom")




## Maps ----
sp = c("Acanthis_flammea","Prunella_collaris") # c score = 0, but low correlation
sp = c("Accipiter_nisus","Buteo_buteo") # c score = 1, rho rather low
dataset_id =5
files = dir('data/occ_5', full.names = T)
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
dat = rbindlist(tmp)
grid = vect(paste0("data/all_scales_atlas_", dataset_id, ".gpkg"))
tmp = dat[dat$scientificName %in% sp,]
p1 = merge(grid, tmp[,.(pres.abs, mean.psi, scientificName, endYear, siteID, scalingID)], all.y=TRUE)
p1 = st_as_sf(p1)
ggplot(p1[p1$scalingID==1,], aes(fill=pres.abs))+
  geom_sf(alpha=0.6)+
  facet_grid(scalingID+scientificName~endYear)
ggplot(p1[p1$scalingID==1,], aes(fill=mean.psi))+
  geom_sf(alpha=0.6)+
  facet_grid(scalingID+scientificName~endYear)


