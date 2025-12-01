# co-occurrence change data per species pair
rm(list=(ls())) # clear workspace
library(groundhog)
pkgs = c("data.table",
         "ggplot2")
groundhog.library(pkgs, "2025-10-25")

old = readRDS("data/processed_spass.rds")
new = readRDS("data/processed_spass2.rds")

nrow(old)
nrow(new)

table(old$dataset_id)
table(new$dataset_id)

table(old$dataset_id, old$overlap)
table(new$dataset_id, new$overlap)

table(old$species_pair %in% new$species_pair)
# that is a quarter, more than 10% hmmmm
sort_pairs_vectorized = function(x) {
  sapply(x, function(item) {
    spl = unlist(strsplit(item, split = "|", fixed = TRUE))
    paste0(sort(spl), collapse = "|")
  })
}

new$species_pair = sort_pairs_vectorized(new$species_pair)
old$species_pair = sort_pairs_vectorized(old$species_pair)
table(old$species_pair %in% new$species_pair)/nrow(old)


table(new$species_pair %in% old$species_pair) 


ggplot(old, aes(x=delta_cor_obs))+
  geom_histogram()+
  facet_wrap('dataset_id', scales="free_y")

ggplot(new, aes(x=delta_cor_obs))+
  geom_histogram()+
  facet_wrap('dataset_id', scales="free_y")


all(new$species_pair %in% old$species_pair)


test <- merge(unique(old[, .(dataset_id, species_pair, delta_cor_obs)]), 
              unique(new[, .(dataset_id, species_pair, delta_cor_obs)]), 
              all.x=TRUE, by=c("dataset_id", "species_pair"))


ggplot(test, aes(x=delta_cor_obs.x, y=delta_cor_obs.y))+
  geom_point()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)

# now check the z-scores
test <- merge(unique(old[time_bin=='T1', .(dataset_id, species_pair, cor_obs, cor_ses)]), 
              unique(new[time_bin=='T1', .(dataset_id, species_pair, cor_obs, cor_ses)]), 
              all.x=TRUE, by=c("dataset_id", "species_pair"))
ggplot(test, aes(x=cor_obs.x, y=cor_obs.y))+
  geom_point(alpha=0.1)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  facet_wrap('dataset_id')
# correlation is exactly the same
ggplot(test, aes(x=cor_ses.x, y=cor_ses.y))+
  geom_point(alpha=0.1)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  facet_wrap('dataset_id')+
  coord_cartesian(xlim=c(-3,+3), ylim=c(-3,+3))
# z-score varies because i did not set a seed, thats ok

fig_hist = test 
fig_hist$x_bin_old = cut(
  fig_hist$cor_ses.x,
  breaks = c(-Inf, -1.96, 1.96, Inf),
  labels = c("Seg (< -1.96)", "Neutral (-1.96 to 1.96)", "Agg (> 1.96)"),
  ordered_result = TRUE   # ensures ordering
)
fig_hist$x_bin_new = cut(
  fig_hist$cor_ses.y,
  breaks = c(-Inf, -1.96, 1.96, Inf),
  labels = c("Seg (< -1.96)", "Neutral (-1.96 to 1.96)", "Agg (> 1.96)"),
  ordered_result = TRUE   # ensures ordering
)

ggplot(fig_hist, aes(x=cor_obs.x))+
  geom_histogram(binwidth = 0.1, aes(fill=x_bin_old), bins=160)+ # aes(y = ..density..), , fill="grey80"
  facet_grid(~dataset_id, scales="free_y")+
  geom_vline(xintercept = 0, lty=1, col="grey40", lwd=0.3)

ggplot(fig_hist, aes(x=cor_obs.y))+
  geom_histogram(binwidth = 0.1, aes(fill=x_bin_new), bins=160)+ # aes(y = ..density..), , fill="grey80"
  facet_grid(~dataset_id, scales="free_y")+
  geom_vline(xintercept = 0, lty=1, col="grey40", lwd=0.3)




told = readRDS("transitions_old.rds")
tnew = readRDS("transitions_new.rds")
test = merge(told,tnew, all=TRUE, by=c("species_pair", "atlas"))

trans1 = as.data.table(table(test$atlas, test$transition.x, useNA = "ifany"))
trans2 = as.data.table(table(test$atlas, test$transition.y, useNA = "ifany"))


table(is.na(test$transition.y))/nrow(test)

ggplot(test, aes(y=transition.x))+
  geom_bar()
ggplot(test, aes(y=transition.y))+
  geom_bar()

# which transitions?
ggplot(test[which(is.na(transition.y)),], aes(y=transition.x))+
  geom_bar()
