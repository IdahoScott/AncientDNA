#PRJNA473816 Baums Acropora palmata & cervicornis data
export BioProject=PRJNA473816
$HOME/edirect/esearch -db sra -query $BioProject | efetch --format runinfo |cut -d "," -f 1 | grep SRR > $BioProject.SRR && $HOME/edirect/esearch -db sra -query $BioProject | efetch --format runinfo > $BioProject.fullMeta.csv

for A in `cat todo`;do
fastq-dump $A; done
ls5_launcher_creator.py -j gets -n gets -a tagmap -e cbscott@utexas.edu -t 04:00:00 -N 8 -w 1 -q normal

ls *.fastq > done
sed 's/.fastq//' done | sort > done.sort
sort unfinished > unfinished.sort
grep -Fxvf done.sort unfinished.sort > todo

export GENOME_FASTA=$WORK/ancientDNAdb/Apalm_and_sym.fasta #super low alignment rate?

>ancientremap #This mapping will have soft clipping!
for F in *.fastq; do
echo "bowtie2 --no-unal --local -x $NEWREF -U $F -S ${F/.fastq}.softclip.sam && \
samtools sort -O bam -o ${F/.fastq}.softclip.sorted.bam ${F/.fastq}.softclip.sam && samtools index -c ${F/.fastq}.softclip.sorted.bam" >> ancientremap; done
ls5_launcher_creator.py -j ancientremap -n ancientremap -a tagmap -e cbscott@utexas.edu -t 03:30:00 -N 1 -w 4 -q normal
sbatch ancientremap.slurm

>fix1766 #This mapping will have soft clipping!
for F in S17466*.fastq; do
echo "bowtie2 --no-unal --local -x $GENOME_FASTA -U $F -S ${F/.fastq}.anc.softclip.sam && \
samtools sort -O bam -o ${F/.fastq}.softclip.sorted.bam ${F/.fastq}.anc.softclip.sam && samtools index -c ${F/.fastq}.softclip.sorted.bam" >> fix1766; done
ls5_launcher_creator.py -j fix1766 -n fix1766 -a tagmap -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 1 -q normal
sbatch fix1766.slurm


>modern_remap
for F in *.fastq; do
echo "bowtie2 --no-unal --local -x $NEWREF -U $F -S ${F/.fastq}.softclip.sam && \
samtools sort -O bam -o ${F/.fastq}.softclip.sorted.bam ${F/.fastq}.softclip.sam && samtools index -c ${F/.fastq}.softclip.sorted.bam" >> modern_remap; done
ls5_launcher_creator.py -j modern_remap -n modern_remap -a tagmap -e cbscott@utexas.edu -t 15:00:00 -N 4 -w 12 -q normal
sbatch modern_remap.slurm

>modern_sort
for F in *.softclip.sam; do
echo "samtools sort -O bam -o ${F/.sam}.sorted.bam $F && samtools index -c ${F/.sam}.sorted.bam" >> modern_sort; done
ls5_launcher_creator.py -j modern_sort -n modern_sort -a tagmap -e cbscott@utexas.edu -t 00:30:00 -N 1 -w 12 -q normal
sbatch modern_sort.slurm


for F in `cat cervicornis`; do
echo "bowtie2 --no-unal --local -x $GENOME_FASTA -U $F -S ${F/.fastq}.default.sam && \
samtools sort -O bam -o ${F/.fastq}.default.sorted.bam ${F/.fastq}.default.sam && samtools index -c ${F/.fastq}.default.sorted.bam" >> maps_modern; done
ls5_launcher_creator.py -j maps_modern -n maps_modern -a mega2014 -e cbscott@utexas.edu -t 20:00:00 -N 2 -w 12 -q normal
sbatch maps_modern.slurm



>sort
for F in *.default.sam; do
echo "samtools sort -O bam -o ${F/.default.sam}.default.sorted.bam $F && samtools index -c ${F/.default.sam}.default.sorted.bam" >> sort; done

ls5_launcher_creator.py -j sort -n sort -a mega2014 -e cbscott@utexas.edu -t 02:00:00 -N 1 -w 12 -q development

ls *.default.sorted.bam* > files
>copy
for file in `cat files`; do
echo "cp $file .." >> copy; done
ls5_launcher_creator.py -j copy -n copy -a mega2014 -e cbscott@utexas.edu -t 00:05:00 -N 1 -w 4 -q normal



#get truncated files
> babyfiles
for file in *.default.sam; do
echo "samtools view -bh $file | head -7000000 > ${file/.default.sam}.default.short.sam" >> babyfiles; done
ls5_launcher_creator.py -j babyfiles -n babyfiles -a tagmap -e cbscott@utexas.edu -t 00:10:00 -N 1 -w 24 -q normal

>babysort
for file in *.default.short.sam; do
echo "samtools sort -O bam -o ${file/.sam}.sorted.bam $file && samtools index -c ${file/.sam}.sorted.bam" >> babysort; done
ls5_launcher_creator.py -j babysort -n babysort -a tagmap -e cbscott@utexas.edu -t 00:20:00 -N 1 -w 24 -q normal

>babycorals
for file in *.default.sorted.bam; do
echo "cat Coral_scaffolds | tr '\n' ' ' | xargs samtools view -bh $file > ${file/.default.sorted.bam}.coral.bam" >> babycorals; done
ls5_launcher_creator.py -j babycorals -n babycorals -a mega2014 -e cbscott@utexas.edu -t 00:10:00 -N 1 -w 12 -q development

#these should have much higher coverage than the *tiny* files (tiny files have about 20,000 reads, short files have about 7,000,000)
cat Coral_scaffolds | tr '\n' ' ' | xargs samtools view -bh S17466.softclip.sorted.bam > S17466.softclip.sorted.coral.bam

#make tiny test files
for file in *.fastq; do
head -80000 $file > ${file}.tiny.fastq; done

>tiny_map
for F in *.tiny.fastq; do
echo "bowtie2 --no-unal --local -x $GENOME_FASTA -U $F -S ${F/.fastq.tiny.fastq}.tinylocal.sam && \
samtools sort -O bam -o ${F/.fastq.tiny.fastq}.tinylocal.sorted.bam ${F/.fastq.tiny.fastq}.tinylocal.sam && samtools index -c ${F/.fastq.tiny.fastq}.tinylocal.sorted.bam" >> tiny_map; done
ls5_launcher_creator.py -j tiny_map -n tiny_map -a tagmap -e cbscott@utexas.edu -t 00:10:00 -N 1 -w 10 -q development


#don't worry about this- something funky is going on with adapters.
cutadapt --format fastq -q 15,15 -m 30 -a AGATCGGA -o SRR7235987.tiny.trim SRR7235987.fastq.tiny.fastq > runlog_cutadapt.txt
cutadapt --format fastq -q 15,15 -m 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o SRR7235987.tinylong.trim SRR7235987.fastq.tiny.fastq > runlog_cutadapt.txt

#tests
bowtie2 --no-unal -x $GENOME_FASTA -U SRR7235987.tinylong.trim -S SRR7235987.tinylong.trim.sam
bowtie2 --no-unal --local -x $GENOME_FASTA -U SRR7235987.tiny.trim -S SRR7235987.tiny.local.trim.sam #fixes the issue!
bowtie2 --no-unal --local -x $GENOME_FASTA -U SRR7235987.fastq.tiny.fastq -S SRR7235987.tiny.local.trim.sam
#try mapping them with --local (see if this returns better results)

#NOW want only coral reads
for file in *.sorted.bam; do
cat Coral_scaffolds | tr "\n" " " | xargs samtools view -bh $file > ${file/.sorted.bam}.coral.bam && echo $file && samtools view -c ${file/.sorted.bam}.coral.bam; done

>coralsub
for file in *.sorted.bam; do
echo "cat Coral_scaffolds | tr '\n' ' ' | xargs samtools view -bh $file > ${file/.sorted.bam}.coral.bam" >> coralsub; done
ls5_launcher_creator.py -j coralsub -n coralsub -a mega2014 -e cbscott@utexas.edu -t 00:30:00 -N 1 -w 12 -q normal


>resort
for file in *.coral.bam; do
echo "samtools sort -O bam -o ${file/.coral.bam}.coral.sorted.bam $file && samtools index -c ${file/.coral.bam}.coral.sorted.bam" >> resort; done
ls5_launcher_creator.py -j resort -n resort -a mega2014 -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 12 -q development




cp /scratch/06909/cbscott/ancientDNA/aDNA_check/map2_coral_zoox/*.default.sorted.coral.bam* .
for file in *.default.coral.bam; do
samtools sort -O bam -o ${file/.default.coral.bam}.ancient.coral.sorted.bam $file && samtools index -c ${file/.default.coral.bam}.ancient.coral.sorted.bam; done


FILTERSQ="-uniqueOnly 1 -remove_bads 1 -minMapQ 20"
TODOQ="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

#try a scaffold name that actually exists, this works
echo "ls *.coral.sorted.bam > bams && angsd -b bams -r 'Sc0a5M3_1;HRSCAF=1': -GL 1 $FILTERSQ $TODOQ -P 12 -out dd && Rscript ~/bin/plotQC.R dd >qualRanks">a0
ls5_launcher_creator.py -j a0 -n a0 -a mega2014 -e cbscott@utexas.edu -t 02:00:00 -w 1 -q development

for file in bams_*; do
  echo "angsd -b $file -r 'Sc0a5M3_1;HRSCAF=1': -GL 1 $FILTERSQ $TODOQ -P 12 -out ${file}.dd && Rscript ~/bin/plotQC.R ${file}.dd > ${file}.qualRanks">>a0;
done
ls5_launcher_creator.py -j a0 -n a0 -a mega2014 -e cbscott@utexas.edu -t 00:10:00 -w 4 -q normal




#need several lists of bams
cp bams bams.qc
#Not going to worry about trimming
export GENOME_REF=$WORK/ancientDNAdb/Apalm_and_sym.fasta
export MinIndPerc=0.5 #picked from dd.pdf - low because bad quality? How to pick this..... Double check this
#we like to keep at at 0.8
#restrict to sites covered in ancient sample ONLY: can do this angsd way... just the bam for the ancient...
#relaxed filters, no other filters, and then extract sites from there. How to extract sites from 'doMaf' allele freq table; cut first four columns
#this becomes 2col table.... angsd sites index that file name... and then use this as -sites argument

#angsd doMath with just one sample... do maj minor, etc. cannibalize the script from 2bRAD pipeline... remove all filters mapping and baseq


JUST_ANC='-minMapQ 20 -minQ 20 -doMaf 1 -doMajorMinor 1 -GL 1 -out anc_test'

#Similar to approach for selecting sites which pass filters.

# initial IBS production, detecting and removing clones (see hctree.pdf and resulting bams.nr)
#bams.qc only has NA values in it? why?
FILTERS0='-minInd $MI -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -doHWE 1 -hetbias_pval 1e-3 -skipTriallelic 1'
TODO0='-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doPost 1 -doGlf 2'
echo 'export NIND=`cat bams.qc | wc -l`; export MI=`echo "($NIND*$MinIndPerc+0.5)/1" | bc`' >calc1
echo "source calc1 && angsd -b bams.qc -GL 1 $FILTERS0 $TODO0 -P 12 -out myresult && Rscript ~/bin/detect_clones.R bams.qc myresult.ibsMat 0.15">a1
ls5_launcher_creator.py -j a1 -n a1 -a mega2014 -e cbscott@utexas.edu -t 02:00:00 -w 1 -q development #might need more time, even for tiny!



echo "angsd -b bams.qc -GL 1 -doMaf 2 -doMajorMinor 1 -P 12 -out nofilt && Rscript ~/bin/detect_clones.R bams.qc nofilt.ibsMat 0.15">a1
echo "angsd -b bams.qc -GL 1 -P 12 $TODO0 -out nofilt && Rscript ~/bin/detect_clones.R bams.qc nofilt.ibsMat 0.15">a1
ls5_launcher_creator.py -j a1 -n a1 -a mega2014 -e cbscott@utexas.edu -t 04:00:00 -w 1 -q normal #might need more time, even for tiny!
#./angsd  -doMaf 2 -doMajorMinor 1 -out TSK -bam bam.filelist -GL 1 -r 1: from website for nofilter angsd

#Now we'll do this for every sample
#forgot the fake clone!
echo S17464.softclip.coral.sorted.bam > S17464
echo S17465.softclip.coral.sorted.bam > S17465
echo S17466.softclip.coral.sorted.bam > S17466
echo S17463.softclip.coral.sorted.bam > S17463

#could make a loop
export file="S17463"
JUST_ANC="-minMapQ 20 -minQ 20 -doMaf 1 -doMajorMinor 1 -GL 1 -out anc_${file}"
angsd -b $file $JUST_ANC
gunzip anc_${file}.mafs.gz
cut -f 1,2,3,4 anc_${file}.mafs | tail -n +2 > sites_${file}.txt #need only first four columns, without a header
cut -f1 sites_${file}.txt | sort | uniq >chrs_${file}.txt
angsd sites index sites_${file}.txt


#Think about also making a chr file... that way we don't analyze regions with no sites!
#Doesn't look great, but let's just try it anyways
#Might need to remap...
export GENOME_REF=$WORK/ancientDNAdb/Apalm_and_sym.fasta
export MinIndPerc=0.5
export file="S17463"
export bams="bams_17463"
export out="myresult17463"
#cp $bams bams.qc
FILTERS0="-minInd $MI -minMapQ 20 -minQ 20 -minMaf 0.05 -dosnpstat 1 -doHWE 1 -hetbias_pval 1e-3 -skipTriallelic 1"
TODO0="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doPost 1 -doGlf 2 -doQsDist 1 -sites sites_${file}.txt -rf chrs_${file}.txt"
echo 'export NIND=`cat $bams | wc -l`; export MI=`echo "($NIND*$MinIndPerc+0.5)/1" | bc`' >calc1
echo "source calc1 && angsd -b $bams -GL 1 $FILTERS0 $TODO0 -P 12 -out $out && Rscript ~/bin/detect_clones.R $bams ${out}.ibsMat 0.15 && cp quality.txt quality_${file}.txt">a_${file}
ls5_launcher_creator.py -j a_${file} -n a_${file} -a mega2014 -e cbscott@utexas.edu -t 01:45:00 -w 1 -q normal #might need more time, even for tiny!
sbatch a_${file}.slurm


########remove minMaf for more variable SNPs !!!
#Differences in ancient sample compared to reference? (allele frequency filtering?)
#extinct variation? But useless for PCA?
#smartPCA
#pseudohaploidization script on github, but bad documentation.
#how many PCS to use to compute bootstrap?
#genome wide heterozygosity 0.003
#0.5 G genome size

#Have to go back to the genotypes to project onto PCA.... CANNOT do from distant matrix

#Get admixuture proportions from the ancient and modern samples


#FIRST STEP- let's bootstrap/jacknife/randomsample scaffolds:

export GENOME_REF=$WORK/ancientDNAdb/Apalm_and_sym.fasta
export MinIndPerc=0.5
export file="S17463"
export bams="bams_17463"
export out="myresult17463"
#cp $bams bams.qc
>index
for file in sampled_sites_S17463*.txt; do
echo "cut -d ' ' -f 1,4,5,6 $file > ${file/.txt}.trim.txt && angsd sites index ${file/.txt}.trim.txt && cut -d ' ' -f1 ${file/.txt}.trim.txt > sampled_chrs_${file#sampled_sites_}" >> index; done
ls5_launcher_creator.py -j index -n index -a tagmap -e cbscott@utexas.edu -t 00:30:00 -w 10 -N 1 -q development

#make chrs
cut -f1 sites_${file}.txt


FILTERS0="-minInd $MI -minMapQ 20 -minQ 20 -minMaf 0.039 -dosnpstat 1 -doHWE 1 -hetbias_pval 1e-3 -skipTriallelic 1"
#TODO0="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doPost 1 -doGlf 2 -doQsDist 1 -sites sites_${file}.txt -rf chrs_${file}.txt"
echo 'export NIND=`cat $bams | wc -l`; export MI=`echo "($NIND*$MinIndPerc+0.5)/1" | bc`' >calc1
for i in {1..100}; do
echo "source calc1 && angsd -b $bams -GL 1 $FILTERS0 -doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doPost 1 -doGlf 2 -doQsDist 1 -sites sampled_sites_${file}_${i}.trim.txt -rf sampled_chrs_${file}_${i}.txt -P 12 -out ${out}.${i} && Rscript ~/bin/detect_clones.R $bams ${out}.${i}.ibsMat 0.15">>a_${file};
done
ls5_launcher_creator.py -j a_${file} -n a_${file} -a tagmap -e cbscott@utexas.edu -t 03:30:00 -w 16 -N 4 -q normal #might need more time, even for tiny!
sbatch a_${file}.slurm #something is wrong here....



#====== make sure there's not worse mapping quality to Acropora millepora....
GENOME_REF=/work/06909/cbscott/sym_references/Amil_symABCD.fasta

>maps_amil
for file in *.fastq; do
echo "bowtie2 --no-unal -x $GENOME_REF -U $file -S ${file/.fastq}.amil.sam && \
samtools sort -O bam -o ${file/.fastq/}.amil.sorted.bam ${file/.fastq/}.amil.sam && samtools index -c ${file/.fastq/}.amil.sorted.bam " >> maps_amil;
done
ls5_launcher_creator.py -j maps_amil -n maps_amil -t 03:30:00 -w 4 -N 1 -a tagmap -e cbscott@utexas.edu -q normal
sbatch maps_amil.slurm

for file in *.sorted.bam; do
cat Amil_scaffs | tr "\n" " " | xargs samtools view -bh $file > ${file/.sorted.bam}.coral.bam && echo $file && samtools view -c ${file/.sorted.bam}.coral.bam; done


for file in *.sorted.bam; do
samtools view -bSq 20 $file > ${file/.sorted.bam}.filtered20.bam && \
samtools sort -O bam -o ${file/.sorted.bam}.filtered20.sorted.bam ${file/.sorted.bam}.filtered20.bam && \
samtools index -c ${file/.sorted.bam}.filtered20.sorted.bam; done

for file in *.filtered20.sorted.bam; do
cat Amil_scaffs | tr "\n" " " | xargs samtools view -bh $file > ${file/.sorted.bam}.coral.bam && echo $file && samtools view -c ${file/.sorted.bam}.coral.bam; done

#===========Unused Qual Check=======
APPLY_ANC= '-minMapQ 20 -minQ 20 -doMaf 1 -doMajorMinor 1 pdoCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doPost 1 -doGlf 2 -GL 1 -sites sites.txt -out sites_test'
angsd -b test $APPLY_ANC


export file="S17465"
export bams="bams_17465"
FILTERSQ="-uniqueOnly 1 -remove_bads 1 -minMapQ 20"
TODOQ="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2 -sites sites_${file}.txt"

#try a scaffold name that actually exists, this works
echo "angsd -b $bams -GL 1 $FILTERSQ $TODOQ -P 12 -out sub.dd && Rscript ~/bin/plotQC.R sub.dd >qualRanks && cp quality.txt quality_${file}.txt">a0
ls5_launcher_creator.py -j a0 -n a0 -a mega2014 -e cbscott@utexas.edu -t 00:02:00 -w 1 -q development
sbatch a0.slurm


#======what about overlapping variable sites from within the ancient samples?


#=======ADMIXTURE==========
for K in `seq 2 5` ;
do
NGSadmix -likes myresult17463.9.beagle.gz -K $K -P 12 -o mydata17463_k${K};
done


#======Recreate Scaffolds with numbering========
>subcorals
for file in *.softclip.sorted.bam; do
echo "cat NumberedCoralScaffolds | tr '\n' ' ' | xargs samtools view -bh $file > ${file/.softclip.sorted.bam}.coral.bam" >> subcorals; done
>sortcorals
for file in *.coral.bam; do
echo "samtools sort -O bam -o ${file/.bam}.sorted.bam $file && samtools index -c ${file/.bam}.sorted.bam" >> sortcorals; done
ls5_launcher_creator.py -j sortcorals -n sortcorals -a tagmap -e cbscott@utexas.edu -t 00:30:00 -N 2 -w 20 -q normal
