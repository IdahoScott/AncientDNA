##ancient DNA IBS
library(ggplot2)
library(vegan)
library(pheatmap)
library(strataG)
library(ggfortify)
file="bams_17465"
mat="myresult17465.ibsMat"
# change the line below to reflect location of PopGenClass directory
meta <- read.csv("Apalm_meta_formatted.csv", header = T, stringsAsFactors = F)
#for snp analysis
submeta <- cbind(meta$Species, meta$SRA.Accession)[1:42,]
write.table(submeta, file = "modern_metadata.txt", quote = F, col.names = F, row.names = F)

#quals <- read.table("quality.txt")  only 248 entries? why? 
#colnames(quals) <- c("Run", "Quality")
#quals$Run <- sub("\\..+", "",quals$Run)
#meta <- merge(meta, quals)
# reading metadata; all we need from this now is "site" - place where mice were caught
bams=scan(file,what="character") #changed from bams.nr
bams=sub("\\..+","",bams)

#add a clone
if (file == "bams_17463")
{
  meta <- rbind(meta, meta[meta$SRA.Accession == "SRR7236019", ])
  meta$SRA.Accession[nrow(meta)] <- "cloneSRR7236019" #fake data- always use the same clone 
}

keep <- data.frame(SRA.Accession=bams)
subbed_meta <- merge(meta, keep)
subbed_meta <- subbed_meta[match(keep$SRA.Accession, subbed_meta$SRA.Accession),]
subbed_meta$Ancient <- ifelse(subbed_meta$age != 0, 'Y', 'N')
#remove <- c("SRR7235987", "SRR7235988")
#nona_meta <- subbed_meta[subbed_meta$SRA.Accession %notin% remove,]
labels <- subbed_meta$SRA.Accession
pop <- subbed_meta$Region
age <- subbed_meta$age
ancient <- subbed_meta$Ancient
species <- subbed_meta$Species

# reading IBS matrix
ibs=as.matrix(read.table(mat))
#names(ibs)=c(names, names)
#rownames(ibs) = labels
#colnames(ibs) = labels
rownames(ibs) = c("oldie", labels[2:length(labels)]) #for procrustes
colnames(ibs) = c("oldie", labels[2:length(labels)]) #for procurustes

ibs = ibs[rownames(ibs) != "cloneSRR7236019",colnames(ibs) != "cloneSRR7236019"]
#ibs = ibs[rownames(ibs) != "oldie", colnames(ibs) != "oldie"] #remove old sample


pheatmap(1-ibs)
hc=hclust(as.dist(ibs),"ave")
plot(hc)
# unconstrained ordination (PCoA)
pp_S17463=capscale(ibs~1)
pp_S17465=capscale(ibs~1)

#procrustes
pro.test <- procrustes(pp_S17463, pp_S17465, symmetric = T) #overlay these guys with ggplot. #infer they are similar by being in the same place 

pro3 <- protest(pp_S17463[, c("MDS1", "MDS2")],
                pp_S17465[, c("Axis.1", "Axis.2")],
                scores = "sites", permutations = how(nperm = 9999))


pro3 <- protest(pp_S17463,
                pp_S17465,
                scores = "sites", permutations = how(nperm = 9999))
pro3

#Do it the ggplot way
ctest <- data.frame(rda1=pro.test$Yrot[,1],
           rda2=pro.test$Yrot[,2],xrda1=pro.test$X[,1],
           xrda2=pro.test$X[,2],anc=ancient, species = species)
ggplot(ctest) +
  geom_point(aes(x=rda1, y=rda2, shape=anc, color=species), size = 3) +
  geom_point(aes(x=xrda1, y=xrda2, shape=anc, color = species), size = 3) +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2, color = species),arrow=arrow(length=unit(0.2,"cm"))) + 
  theme_minimal() + ggtitle("Procrustes Test Comparison of 1000 y.o. Ancient Samples")


# simple plot:
plot(pp)  #need the coloration

# eigenvalues: how many good eigenvalues (exceeding broken stick expectation) are there:
plot(pp$CA$eig/sum(pp$CA$eig),ylab="proportion of variance explained",xlab="eigenvector #")
lines(bstick(length(pp$CA$eig)),col="red") # "broken stick" model: expected random distribution of eigenvalues
# looks like we have 2 good eigenvalues

# extracting "scores" table, to plot
axes2plot=c(1,2) # which PCAs to plot
scores=data.frame(pp_S17463$CA$u[,axes2plot])
scores<-scores[labels,]

#quartz()
ggplot(scores,aes(scores[,1],scores[,2],color=species)) + geom_point(alpha=0.8, size = 5)+theme_bw()+xlab(names(scores)[1])+ylab(names(scores)[2])+coord_equal() +
  ggtitle("S17463 PCA")


#Project ancient sample onto this PCA...?
# perform principal components analysis
pca <- prcomp(ibs) 


old=as.matrix(read.table(mat))
rownames(old) = c("oldie", labels[2:length(labels)]) #for procrustes
colnames(old) = c("oldie", labels[2:length(labels)]) #for procurustes
newdata <- old[rownames(old) == "oldie", colnames(old == "oldie")]
newdata <- newdata[names(newdata) != "cloneSRR7236019"]
newdata <- newdata[names(newdata) != "oldie"]
newdata <- t(as.data.frame(newdata))



# project new data onto the PCA space... and add to PCA plot?
S17463 <- scale(newdata, pca$center, pca$scale) %*% pca$rotation
rownames(S17463) <- c("S17463")
add_new <- rbind(pca$x, S17463)

newmeta <- subbed_meta[subbed_meta$SRA.Accession %in% rownames(add_new),]
newmeta <- newmeta[match(keep$SRA.Accession, subbed_meta$SRA.Accession),]
species <- na.omit(newmeta$Species)
ancient <- na.omit(newmeta$Ancient)

#plot(pca$x[,1], pca$x[,2])
ggplot(add_new, aes(add_new[,1], add_new[,2], color = species, shape = ancient)) + geom_point() # looks pretty good. Now add a pont 

#autoplot(pca, data = ibs, color = 'species')

scaling <- pca$sdev[1:2] * sqrt(nrow(ibs))

pc1 <- rowSums(t(t(sweep(ibs[,2:7], 2 ,colMeans(ibs[,2:7]))) * s.eigen$vectors[,1] * -1) / scaling[1])
#pc2 <- rowSums(t(t(sweep(pilots[,2:7], 2, colMeans(pilots[,2:7]))) * s.eigen$vectors[,2]) / scaling[2])
#Collect the PCs in a data.frame and plot using ggplot (loaded when ggfortify was loaded).



#can we accomplish this with vegan?




pp1=capscale(ibs~ancient)
axes2plot=c(1,2)
scores=data.frame(pp1$CA$u[,axes2plot])
ggplot() + geom_point(data=scores,aes(scores[,1],scores[,2],color=pop, shape=ancient),alpha=0.5, size = 3)+theme_bw()+xlab(names(scores)[1])+ylab(names(scores)[2])
#hclust(pp1$CCA$Xbar)
#filter by quality scores?

quals <- read.table("quality.txt") #something is wrong here. qual should not be zero for big files?
quals$V1 <- sub("\\..+","",quals$V1)
colnames(quals) <- c("SRA.Accession", "quality")
qual_meta <- merge(subbed_meta, quals)


#ibs_qual <- ibs_nona[qual_meta$SRA.Accession, qual_meta$SRA.Accession]
q_prob <- qual_meta$quality
qual_ancient <- qual_meta$Ancient
qual_bp <- 2020-qual_meta$age
# unconstrained ordination (PCoA)
pp_q=capscale(ibs~q_prob+qual_ancient)
# simple plot:
plot(pp_q)  #need the coloration

pp1=capscale(ibs_qual~q_prob+qual_ancient)
axes2plot=c(1,2,3)
scores=data.frame(pp1$CA$u[,axes2plot])
ggplot() + geom_point(data=scores,aes(scores[,1],scores[,2],color=qual_ancient),alpha=0.8)+theme_bw()+xlab(names(scores)[1])+ylab(names(scores)[2]) +ggtitle("Constrained PCA: ibs~quality+ancient")
library(plotly)
plot_ly(x= scores[,1], y=scores[,2], z=scores[,3], color = qual_ancient)
test <- dist(pp1$CA$u)
hc = hclust(test) #is this eve?n legal???? Stats???? Seems like fair game. PCA give distance in eigen space, can we map that to cluster dendrogram....

library(factoextra)
library(FactoMineR)

res.pca <- PCA(pp1$CA$u, ncp = 5, graph = TRUE)
res.hcpc <- HCPC(res.pca, graph = T)
fviz_dend(res.hcpc, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8      # Augment the room for labels
)
plot(res.hcpc, choice = "3D.map") #ridiculogram of the day

