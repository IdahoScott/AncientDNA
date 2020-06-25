#Distribution of PMD scores for reads from one sample

#At some point want to know PMD scores for reads from coral and syms only

pmdscores <- read.table("PMD_dist.txt") #the file has 4 columns... not sure what the first three are
hist(pmdscores[,4])

P = ecdf(pmdscores[,4])
plot(P, xlim=c(-3, 8), ylab= "Percent Reads NOT kept (lines at 90% and 95%)", xlab = "PMD Score, 6000 y/o human aDNA cutoff=3")

abline(v=0, lty = 'dashed', col = 'red')
abline(v=3, lty = 'dashed', col = 'red')
abline(v=1.5, lty = 'dashed', col ='red')
abline(h = 0.9, col = 'blue')
abline(h = 0.95, col = 'blue')
