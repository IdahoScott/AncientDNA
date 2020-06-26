#Distribution of PMD scores for reads from one sample

#At some point want to know PMD scores for reads from coral and syms only
sample_scores <- list.files(pattern = "\\.sortedpmdscores.txt")
pmdscores <- list()
for (i in 1:length(sample_scores)){
pmdscores[[i]] <- read.table(sample_scores[i]) #the file has 4 columns... not sure what the first three are
}
hist(pmdscores[[1]][,4])

P = ecdf(pmdscores[[1]][,4])
P1 = ecdf(pmdscores[[2]][,4])
P2 = ecdf(pmdscores[[3]][,4])
P3 = ecdf(pmdscores[[4]][,4])
plot(P, xlim=c(-3, 8), lwd = 0.5, ylab= "Percent Reads NOT kept (lines at 90% and 95%)", xlab = "PMD Score, 6000 y/o human aDNA cutoff=3", main = "Distribution of PMD scores")
plot(P1, add = T, col = 'blue4', lwd = 0.5)
plot(P2, add = T, col = 'red4', lwd = 0.5)
plot(P3, add = T, col = 'green4', lwd = 0.5)
abline(v=0, lty = 'dashed', col = 'red')
abline(v=3, lty = 'dashed', col = 'red')
abline(v=1.5, lty = 'dashed', col ='red')
abline(h = 0.9, col = 'blue')
abline(h = 0.95, col = 'blue')
