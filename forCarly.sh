##First install some basic dependencies
#Load the intel compiler, gsl and python3 modules
module load intel/18.0.2
#module load gsl/2.3
module load gsl
module load python3/3.7.0
module load samtools
module load Rstats/4.0.3

#Install anaconda in a folder that has enough space (4-5Gb).
wget https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh #put in work
bash Anaconda3-2020.02-Linux-x86_64.sh #[ACCEPT THE LICENSE AGREEMENT, INSTALL THIS IN A FOLDER AND THEN INITIALIZE IT WHEN PROMPTED]
#Install openblas
$WORK/anaconda3/bin/conda install -c anaconda openblas #[THE LOCATION OF WHERE YOU INSTALLED CONDA WILL BE DIFFERENT FROM MINE]


#I think this is the pseudohaploidization
#Run pileupCaller. Please see (https://github.com/stschiff/sequenceTools) for the various arguments and input:
#This program should output 3 files, in Eigenstrat format. See more here: https://reich.hms.harvard.edu/software/InputFileFormats for the file formats
samtools mpileup -R -B -q10 -Q20 -l /work/07475/vagheesh/stampede2/software/aDnaTools/1240kSNP.bed -f /work/07475/vagheesh/stampede2/software/aDnaTools/hs37d5.fa /work/07475/vagheesh/stampede2/forOthers/forCarly/I6113.merged.d.bam | /work/07475/vagheesh/stampede2/software/aDnaTools/pileupCaller --randomHaploid --sampleNames I6113 --samplePopName Harappan -f /work/07475/vagheesh/stampede2/software/aDnaTools/1240kSNP.snp -e /work/07475/vagheesh/stampede2/forOthers/forCarly/harappan

# list_of_positions /work/07475/vagheesh/stampede2/software/aDnaTools/1240kSNP.bed
#what is a bed file?? #need scaffold name, scaffold start, and scaffold end
#Use bamtobed from bedtools?
#going to have to rename and remap.
awk '/^>/{print ">" ++i; next}{print}' < renamedApalmata.fasta > numberedApalmata.fasta
cat $STOCKYARD/sym_references/symABCD_longref.fasta numberedApalmata.fasta > numberedApalm_and_sym.fasta
export NEWREF=/work/06909/cbscott/lonestar/ancientDNAdb/mastergenome/numberedApalm_and_sym.fasta

echo "bowtie2-build $NEWREF $NEWREF" >btb
ls5_launcher_creator.py -j btb -n btb -l btbl -t 02:00:00 -a tagmap -e cbscott@utexas.edu -w 1 -q normal
sbatch btbl
samtools faidx $NEWREF

#map with soft clipping & then only keep coral reads #see compareWGS script for more details

#break it down line by line
samtools mpileup -R -B -q10 -Q20 -l /work/07475/vagheesh/stampede2/software/aDnaTools/1240kSNP.bed #scaffoldname and positions (should have this output from ANGSD?)
-f /work/07475/vagheesh/stampede2/software/aDnaTools/hs37d5.fa #this is the reference genome
/work/07475/vagheesh/stampede2/forOthers/forCarly/I6113.merged.d.bam #just my aligned bam |
/work/07475/vagheesh/stampede2/software/aDnaTools/pileupCaller --randomHaploid --sampleNames I6113 #name of my sample (should come from the name of my bam)
 --samplePopName Harappan
-f /work/07475/vagheesh/stampede2/software/aDnaTools/1240kSNP.snp -e /work/07475/vagheesh/stampede2/forOthers/forCarly/harappan

#let's see what the outputs are from just samtools mpileup
#should the first argument just be restricted to sites in the ancient file?
export GENOME_REF=/work/06909/cbscott/lonestar/ancientDNAdb/mastergenome/numberedApalm_and_sym.fasta
#export SITES=/scratch/06909/cbscott/ancientDNA/aDNA_check/compareWGS/all_coral_bams/sites_S17463_reduced.txt
export BED=/scratch/06909/cbscott/ancientDNA/aDNA_check/VagheeshPipeline/CorrectlyNumbered/numberedApalm.bed #if we create this from bams, it will actually limit to our regions of interest
export SITES=/scratch/06909/cbscott/ancientDNA/aDNA_check/VagheeshPipeline/CorrectlyNumbered/vcftest/allanc.snp
export BAM=/scratch/06909/cbscott/ancientDNA/aDNA_check/VagheeshPipeline/CorrectlyNumbered/S17463.softclip.sorted.bam
export bam_list=/scratch/06909/cbscott/ancientDNA/aDNA_check/VagheeshPipeline/CorrectlyNumbered/bamlist
export modernbams=/scratch/06909/cbscott/ancientDNA/aDNA_check/VagheeshPipeline/CorrectlyNumbered/modernbams
export outpref=$SCRATCH/ancientDNA/aDNA_check/VagheeshPipeline/CorrectlyNumbered/anc17463test
export SNP=$SCRATCH/ancientDNA/aDNA_check/VagheeshPipeline/CorrectlyNumbered/vcftest/allanc.snp
#export BED=$SCRATCH/ancientDNA/aDNA_check/VagheeshPipeline/Apalm_and_sym.bed
samtools mpileup -R -B -q10 -Q20 -l $BED -f $GENOME_REF $BAM | /work/07475/vagheesh/stampede2/software/aDnaTools/pileupCaller --randomHaploid --sampleNames Anc17463 -f $SNP -e $outpref
#samtools mpileup -R -B -q10 -Q20 -l $SITES -f $GENOME_REF -b $bam_list

#I ALSO NEED TO DO THIS TO THE MODERN FILES IN ORDER TO CREATE THE INPUT I WILL NEED LATER FOR THE PCA


#potentially need to replace all of my scaffold names with numbers (this might be helpful broadly?)
#try to create .snp file
echo "bcftools mpileup -f $NEWREF -b $bam_list --threads 10 | bcftools call -mv -Ob -o calls.bcf" > variantjob #assumes diploid. this is okay.
ls5_launcher_creator.py -j variantjob -n variantjob -a tagmap -e cbscott@utexas.edu -t 01:00:00 -N 1 -q development

#gotta parallelize bcf tools
#Create subsets of the scaffolds
line1=100
line2=200
tail -n +$line1 NumberedCoralScaffolds | head -n $((line2-line1)) > second_scaffs
line3=200
line4=300
tail -n +$line3 NumberedCoralScaffolds| head -n $((line4-line3)) > third_scaffs
head -99 NumberedCoralScaffolds > first_scaffs
tail -142 NumberedCoralScaffolds > fourth_scaffs
>paralellized_calls
echo "parallel -a first_scaffs 'bcftools mpileup -Ou -f $GENOME_REF -r {} -b $modernbams | bcftools call -m -v -Oz -o {}.vcf.gz'" > parallelized_calls
echo "parallel -a second_scaffs 'bcftools mpileup -Ou -f $GENOME_REF -r {} -b $modernbams | bcftools call -m -v -Oz -o {}.vcf.gz'" >> parallelized_calls
echo "parallel -a third_scaffs 'bcftools mpileup -Ou -f $GENOME_REF -r {} -b $modernbams | bcftools call -m -v -Oz -o {}.vcf.gz'" >> parallelized_calls
echo "parallel -a fourth_scaffs 'bcftools mpileup -Ou -f $GENOME_REF -r {} -b $modernbams | bcftools call -m -v -Oz -o {}.vcf.gz'" >> parallelized_calls
ls5_launcher_creator.py -j parallelized_calls -n parallelized_calls -a tagmap -e cbscott@utexas.edu -t 03:00:00 -N 4 -w 1 -q normal
#still running out of memory for block 441.... break into 4 chunks instead of 2?
#be sure to gunzip things!
>properorder
for i in {1..441}; do
echo "${i}.vcf" >> properorder
done

#combine files
cat properorder | tr '\n' ' ' | xargs bcftools concat > allsamples_anc_modern.vcf
#now i think we need to get rid of anything but biallelic snps
bcftools view -m2 -M2 -v snps allsamples_anc_modern.vcf > allsamples_anc_modern_snps.vcf
#create .snp file
sed '/^#/d' allsamples_anc_modern_snps.vcf > noheader_allsamples_anc_modern.txt
cut -f1,2,4,5 noheader_allsamples_anc_modern.txt > requiredcolsvcf.txt
wc -l requiredcolsvcf.txt #need for next step  1320763
#create column of ids
>snpids
>genpos
for i in {1..1180285}; #this is not short.... should probably be a job maybe faster with awk?
do
    echo "allanc${i}" >> snpids;
    echo "0" >> genpos;
done

cut -f1 requiredcolsvcf.txt > snpchr
cut -f2,3,4 requiredcolsvcf.txt > snpposrefalt
#finally! We have the SNP file - I think...
paste snpids snpchr genpos snpposrefalt > allanc.snp

#create a bed file
#We have this in our reference genome index!
cut -f1,2,3 numberedApalm_and_sym.fai > numberedApalm_and_sym.bed

#To create an .snp formatted file:
#col 1 - SNP name
#col 2 - chromosome
#col 3 - genetic position in morgans, if known. If not known, set to 0.0
#col 4 - physical position in bases
#optional col 5, col 6 for reference and variant alleles

###Now I think I just want to run pilupCaller on one file? on all files?

#NEED TO edit metadata to create these genotype files for modern samples.
#Created the file in R....

#SOMETHING IS WRONG WITH THE METADATA INFO!!!
cut -d' ' -f2 modern_metadata.txt > modernpops
cut -d' ' -f3 modern_metadata.txt > modernnames
awk 'NF{print $0 ".coral.sorted.bam"}' modernnames > modernbams #accidentally overwrote file- should be fine it's just reordered

echo "samtools mpileup -R -B -q10 -Q20 -l $BED -f $GENOME_REF -b $modernbams | /work/07475/vagheesh/stampede2/software/aDnaTools/pileupCaller --randomHaploid --sampleNames `cat modernnames | tr '\n' ','` -f $SNP " > modtest #don't use -e modern option.... what does this create?
sed 's/.coral.sorted.bam$//' modernbams > correct_modnames
echo "samtools mpileup -R -B -q10 -Q20 -l $BED -f $GENOME_REF -b $modernbams | /work/07475/vagheesh/stampede2/software/aDnaTools/pileupCaller --randomHaploid --sampleNames `cat correct_modnames | tr '\n' ',' | sed 's/.$//'` -f $SNP -e modern" > modtest #don't use -e modern option.... what does this create?

#be sure to edit it to get rid of extra comma at the end of the lists.... could probably edit the input files instead?
#hm giving it the pop gave funky results.... remove and restart. Might have to run separately for different pops then merge?

ls5_launcher_creator.py -j modtest -n modtest -a tagmap -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 1 -q development


##PCA Analysis
#Create a new geno file, this time with 22 chromsomes with one of them deleted in turn
a <- read.table("/work/07475/vagheesh/stampede2/forOthers/forCarly/harappan.snp.txt",sep="\t", header=F)
b <- read.fwf("/work/07475/vagheesh/stampede2/forOthers/forCarly/harappan.geno.txt",widths=c(1))
c <- cbind(a,b)
for (j in 1:22)
{
	c <- cbind(c,b)
}
for (j in 1:22)
{
	c[which(c[,2]==j),8+(j-1)] <- 9
}
c$x <- apply(c[ ,c(7:29)],1,paste, collapse = "")
d <- c$x
write.table(d, "/work/07475/vagheesh/stampede2/forOthers/forCarly/harappan.missing.geno", sep="\t", row.names=F, col.names=F, quote=F)
#Create a new ind file reflecting changes to each chromosome
a <- read.table("/work/07475/vagheesh/stampede2/forOthers/forCarly/harappan.ind.txt", sep="\t", header=F)
a$V1 <- as.character(a$V1)
a$V2 <- as.character(a$V2)
a$V3 <- as.character(a$V3)
for (j in 1:22)
{
	a <- rbind(a, c(paste(a$V1[1],".",j,sep=""), "U", paste(a$V3[1],".",j,sep="")))
}
write.table(a, "/work/07475/vagheesh/stampede2/forOthers/forCarly/harappan.missing.ind", sep="\t", row.names=F, col.names=F, quote=F)
#Finally relabel the snp file
cp /work/07475/vagheesh/stampede2/forOthers/forCarly/harappan.snp.txt /work/07475/vagheesh/stampede2/forOthers/forCarly/harappan.missing.snp
#Now we merge this dataset with that of the other samples from modern data
/work/07475/vagheesh/stampede2/software/AdmixTools/bin/mergeit -p /work/07475/vagheesh/stampede2/forOthers/forCarly/parMergeIt #I ran into errors....
mergeit -p parMergeIttest #works
#First start with PCAs of all the samples in the dataset
sbatch /work/07475/vagheesh/stampede2/forOthers/forCarly/pcaEastEurasia.sh
sed 1d /work/07475/vagheesh/stampede2/forOthers/forCarly/v33.1.HO.harappan.evec | tr -s ' ' | sed -e 's/^[ \t]*//' | awk 'NF{NF-=1};1' | awk -v OFS="\t" '$1=$1' | sed -e '1s/^/Sample_ID\tpc1\tpc2\tpc3\tpc4\tpc5\tpc6\tpc7\tpc8\tpc9\tpc10\n/' > /work/07475/vagheesh/stampede2/forOthers/forCarly/v33.1.HO.harappan.evec.txt
#Now obtain a set of PCA coordinates for just the leave one out samples
grep I6113 /work/07475/vagheesh/stampede2/forOthers/forCarly/v33.1.HO.harappan.evec.txt | tr -s ' ' | sed -e 's/^[ \t]*//' | cut -f1,2,3 | tail -22 > /work/07475/vagheesh/stampede2/forOthers/forCarly/I6113.pca.txt

#Compute the block jacknife
a<-read.table("/work/07475/vagheesh/stampede2/forOthers/forCarly/I6113.pca.txt",sep="\t",header=F)
jackknife <- function(x,print=T)
{
	c <- mean(x)
	d <- 0
	for (j in 1:length(x))
	{
		d<-d + (x[j]-c)*(x[j]-c)
	}
	sqrt(d*21/22)
}
jackknife(a$V2)
jackknife(a$V3)
