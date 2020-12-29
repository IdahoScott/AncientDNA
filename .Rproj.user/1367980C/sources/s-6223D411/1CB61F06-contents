#resampled IBS
library(vegan)

file="bams_17465"
path="S17465"
ibsnames="myresult17465"
age = "1000"
bams=scan(file,what="character") #changed from bams.nr, pull from TACC
bams=sub("\\..+","",bams)
meta <- read.csv("Apalm_meta_formatted.csv", header = T, stringsAsFactors = F) #assoc. metadata. I have this locally- it's an edited version of the paper's metadata
#if (file == "bams_17463") #for this sample, I have an artifical clone in the data as a sanity check
#{
 # meta <- rbind(meta, meta[meta$SRA.Accession == "SRR7236019", ])
 # meta$SRA.Accession[nrow(meta)] <- "cloneSRR7236019" #fake data- always use the same clone 
#}
keep <- data.frame(SRA.Accession=bams)
subbed_meta <- merge(meta, keep)
subbed_meta <- subbed_meta[match(keep$SRA.Accession, subbed_meta$SRA.Accession),] #merge sometimes reorders samples, this puts them back properly
subbed_meta$Ancient <- ifelse(subbed_meta$age != 0, 'Y', 'N')
ancient <- subbed_meta$Ancient
species <- subbed_meta$Species

# read in all new IBS matrices 
resample_list <- list()
for(i in 1:100) { #where 100 is the number of random resamplings
  resample_list[[i]] <- read.table(paste("resampled_IBS/", path, "/", ibsnames, ".", i, ".ibsMat", sep = ''), stringsAsFactors = F)
  colnames(resample_list[[i]]) <- bams
  rownames(resample_list[[i]]) <- bams
}
#some of the resamplings only returned NAs (probably due to lack of variable sites)
#this creates big downstream issues, remove the NAs here
resample_list <- Filter(function(a) any(!is.na(a)), resample_list)

perc_kept <- length(resample_list)/100
print(paste("Percent of resamplings kept:", perc_kept*100, "%"))


#first find pairwise procrustes. Need a list of IBS mats.... 
##Create capscale objects for each ibs matrix
ibs_list <- list()
for (i in 1:length(resample_list)){
  ibs_list[[i]] =capscale(resample_list[[i]]~1)
}

perc_new <- length(ibs_list)/length(resample_list)
print(paste("Percent of resamplings kept at this stage:", perc_new*100, "%"))

#calculate procrustes sum of squares error for each PCA
pro_error <- matrix(NA, nrow = length(ibs_list), ncol = length(ibs_list))
for(i in 1:length(ibs_list)){
  for (j in 1:length(ibs_list)){
    pro_error[i, j] <- procrustes(ibs_list[[i]], ibs_list[[j]], symmetric = T)$ss
  }
}

#reference sample: chosen as the sample with the lowest overall error (or distance) from all other samples
ref_samp <- which.min(colSums(pro_error)) # this is the wrong way to think about it! We want the sample closest to the center..... 

#pro_test list
pro_list <- list()
p_vals <- list()
for (i in 1:length(ibs_list)){
pro_list[[i]] <- procrustes(ibs_list[[ref_samp]], ibs_list[[i]], choices=c(1:10), symmetric = T)
#p_vals[[i]] <- protest(ibs_list[[ref_samp]],
#                       ibs_list[[i]],
 #                       scores = "sites", permutations = how(nperm = 999)) #slowish
}

ctest <- data.frame(rda1=pro_list[[1]]$Yrot[,1],
                    rda2=pro_list[[1]]$Yrot[,2],xrda1=pro_list[[1]]$X[,1],
                    xrda2=pro_list[[1]]$X[,2], strap=rep(1), ss = pro_list[[1]]$ss, anc = ancient, species = species, name = bams)
for (i in 2:length(pro_list)){
  foo <- data.frame(rda1=pro_list[[i]]$Yrot[,1],
                    rda2=pro_list[[i]]$Yrot[,2],xrda1=pro_list[[i]]$X[,1],
                    xrda2=pro_list[[i]]$X[,2], strap=rep(i), ss = pro_list[[i]]$ss, anc = ancient, species = species, name = bams)
  ctest <- rbind(ctest, foo)
}

ggplot(ctest) +
  geom_point(aes(x=rda1, y=rda2, col = species), size = 3, alpha = 0.2) +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2, col = species),arrow=arrow(length=unit(0.2,"cm")), alpha = 0.2) + 
  geom_point(aes(x=xrda1, y=xrda2), size = 3, pch= 1) +
  theme_minimal() + ggtitle("Procrustes Test 1000 y.o. Ancient Samples")

ggplot(ctest) +
  stat_ellipse(aes(rda1, rda2, group = name, fill = interaction(anc,species)),type = "norm", geom = "polygon", alpha = 0.4, level = 0.95) +
  #stat_ellipse(aes(rda1, rda2, group = name),type = "norm", linetype = 1) + 
  #geom_point(aes(x=rda1, y=rda2), size = 1.5, alpha =.5) +
  theme_minimal() + ggtitle(paste("Procrustes Test", age, "y.o. Ancient Sample,", path))

# can we make an averaged procrustes? YIELDS SAME QUALITATIVE RESULT
###THIS CODE ONLY DIFFERS IN THAT A REFERENCE SAMPLE IS NOT CHOSEN
#Instead 
avg_sample <- as.matrix(Reduce("+", resample_list) / length(resample_list))
#just replace NaNs with 1s.... meaning there is NO similarity between NA values
avg_sample[is.nan(avg_sample)] <- 1
avg_ibs <- capscale(avg_sample~1)

pro_list <- list()
for (i in 1:length(ibs_list)){
  pro_list[[i]] <- procrustes(avg_ibs, ibs_list[[i]], symmetric = T)
}

ctest <- data.frame(rda1=pro_list[[1]]$Yrot[,1],
                    rda2=pro_list[[1]]$Yrot[,2],xrda1=pro_list[[1]]$X[,1],
                    xrda2=pro_list[[1]]$X[,2], strap=rep(1), ss = pro_list[[1]]$ss, anc = ancient, species = species, name = bams)
for (i in 2:length(pro_list)){
  foo <- data.frame(rda1=pro_list[[i]]$Yrot[,1],
                    rda2=pro_list[[i]]$Yrot[,2],xrda1=pro_list[[i]]$X[,1],
                    xrda2=pro_list[[i]]$X[,2], strap=rep(i), ss = pro_list[[i]]$ss, anc = ancient, species = species, name = bams)
  ctest <- rbind(ctest, foo)
}

ggplot(ctest) +
  geom_point(aes(x=rda1, y=rda2, col = species), size = 3, alpha = 0.2) +
  #geom_point(aes(x=xrda1, y=xrda2), size = 3) +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2, col = species),arrow=arrow(length=unit(0.2,"cm")), alpha = 0.2) + 
  theme_minimal() + ggtitle("Procrustes Test 1000 y.o. Ancient Samples") 

#TRY ORDIELLIPSE

##Try to unstack samples and then color by sample 
ggplot(ctest) +
  geom_point(aes(x=rda1, y=rda2, col = rownames(ctest)), size = 3, alpha = 0.2) +
  #geom_point(aes(x=xrda1, y=xrda2), size = 3) +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2, col = species),arrow=arrow(length=unit(0.2,"cm")), alpha = 0.2) + 
  theme_minimal() + ggtitle("Procrustes Test 1000 y.o. Ancient Samples") 

#THIS CODE PLAYS WITH ADDING CONTOUR LINES. Ultimately, this is kind of useless. 
ggplot(ctest) +
  geom_point(aes(x=rda1, y=rda2, col = species), size = 3, alpha = 0.2) +
  geom_point(aes(x=xrda1, y=xrda2), size = 3) +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2, col = species),arrow=arrow(length=unit(0.2,"cm")), alpha = 0.2) + 
  theme_minimal() + ggtitle("Procrustes Test 1000 y.o. Ancient Samples") + geom_density_2d(aes(x=rda1, y=rda2), col = "black") 

ggplot(ctest) +
  geom_point(aes(x=rda1, y=rda2, col = species), size = 3, alpha = 0.2) +
  geom_point(aes(x=xrda1, y=xrda2), size = 3) +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2, col = species),arrow=arrow(length=unit(0.2,"cm")), alpha = 0.2) + 
  theme_minimal() + ggtitle("Procrustes Test 1000 y.o. Ancient Samples") + geom_density_2d(aes(x=xrda1, y=xrda2))

