#Admixture for 17465

dir="troubleshoot/" # path to input files
inName="mydata17463_k2.qopt" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run
npops=2

file="bams_17463"
bams=scan(paste(dir,file,sep=''),what="character") #changed from bams.nr, pull from TACC
bams=sub("\\..+","",bams)
meta <- read.csv("Apalm_meta_formatted.csv", header = T, stringsAsFactors = F)
keep <- data.frame(SRA.Accession=bams)
subbed_meta <- merge(meta, keep)
subbed_meta <- subbed_meta[match(keep$SRA.Accession, subbed_meta$SRA.Accession),] #merge sometimes reorders samples, this puts them back properly
subbed_meta$Ancient <- ifelse(subbed_meta$age != 0, 'Y', 'N')
subbed_meta$Species[1] <- "Ancient"
subbed_meta$Region[1] <- "Ancient"

pops= data.frame(cbind(subbed_meta$SRA.Accession, subbed_meta$Species)) # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
tbl=read.table(paste(dir,inName,sep=""),header=F)
#i2p=read.table(paste(dir,pops,sep=""),header=F)
names(pops)=c("ind","pop")
#tbl = head(tbl, 38)
tbl=cbind(tbl,pops)
tbl=head(tbl,38)
row.names(tbl)=tbl$ind

head(tbl,20) # this is how the resulting dataset must look

source("plot_Admixture_v4_function.R")

ords=plotAdmixture(data=tbl[order(tbl$pop),],npops=npops,angle=0,vshift=0,hshift=0)


