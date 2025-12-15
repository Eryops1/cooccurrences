# Script: 05_figure5.R
# Author: Melanie Tietje
# Email: tietje@fzp.czu.cz
# GitHub: @Eryops1
# Purpose: Do the math for schematic figure 5.
# Notes: groundhog will ensure the exact same R packages will be used and create
#        a library on first run, which might take a while.
#        Run this script for each atlas separately by adjusting the dataset_id variable


# joint movement ----------------------------------------------------------

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







# expansion (or contraction) -----------------------------------------------

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









