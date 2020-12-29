#aDNA Check/Working Pipeline
#Background: prior aDNA work yielded different results than this time. We need to
#check that our results are actually valid. This will involve:
#1. Remapping the original fastq files. Potentially recreating the original fastq files using samtools
#1a. What chromosomes are the files mapping to?
#2. Run mapDamage on the output, with unfiltered quality.
#3. Run mapDamage on filtered quality output. Want Q=20 and Q=30.

export MY_DIR=/scratch/06909/cbscott/ancientDNA/aDNA_check

#GET FILES
wget https://www.dropbox.com/sh/vl2vwz1yc9kyc9u/AAAVRQ9RJgAJhzyXQ2IUuxtXa?dl=0 #rename file as a .zip, then extract
cp AAAVRQ9RJgAJhzyXQ2IUuxtXa?dl=0 fastqs.zip
unzip fastqs.zip

#convert back to fastq for remapping
cd 201908_shotgun/
>convert
for file in *.bam; do
echo "samtools fastq $file > ${file/.Y1.E1.L1.bam}.fastq" >> convert; done
ls5_launcher_creator.py -j convert -n convert -t 00:30:00 -w 4 -N 1 -a tagmap -e cbscott@utexas.edu -q development
sbatch convert.slurm
rm convert* #keep this directory cleaner

cd ../map2_coral_zoox
export GENOME_FASTA=$WORK/ancientDNAdb/Apalm_and_sym.fasta #have indexed ref genome here. Be careful, it seems we accidentally mapped to the WRONG genome off the bat.
>maps_default
for file in ../201908_shotgun/*.fastq; do
echo "bowtie2 --no-unal -x $GENOME_FASTA -U $file -S ${file/.fastq}.default.sam && \
samtools sort -O bam -o ${file/.fastq/}.default.sorted.bam ${file/.fastq/}.default.sam && samtools index -c ${file/.fastq/}.default.sorted.bam " >> maps_default;
done
ls5_launcher_creator.py -j maps_default -n maps_default -t 00:45:00 -w 4 -N 1 -a tagmap -e cbscott@utexas.edu -q normal
sbatch maps_default.slurm

>maps_clip
for file in S17465.fastq; do
echo "bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x $GENOME_FASTA -U $file -S ${file/.fastq}.softclip.sam && \
samtools sort -O bam -o ${file/.fastq/}.softclip.sorted.bam ${file/.fastq/}.softclip.sam && samtools index -c ${file/.fastq/}.softclip.sorted.bam " >> maps_clip;
done
ls5_launcher_creator.py -j maps_clip -n maps_clip -t 03:30:00 -w 1 -N 1 -a mega2014 -e cbscott@utexas.edu -q normal
sbatch maps_clip.slurm




#Now let's get perc by type....
#get coral scaffolds from source;
grep ">" $GENOME_FASTA | sed 's/^.//' > all_scaffolds
wc -l all_scaffolds 445
head -441 all_scaffolds > Coral_scaffolds

#Filter every bam for coral reads
for file in *default.sorted.bam; do
cat Coral_scaffolds | tr "\n" " " | xargs samtools view -bh $file > ${file/.sorted.bam}.coral.bam && echo $file && samtools view -c ${file/.sorted.bam}.coral.bam; done

#manually did this for cladeC and cladeA... there is something wrong with zooxType
for file in *default.sorted.bam; do
samtools view -bh $file chr13 > ${file/.sorted.bam}.cladec.bam && echo $file && samtools view -c ${file/.sorted.bam}.cladec.bam; done

#filter by first qual: 20
for file in *.sorted.bam; do
echo " samtools view -bSq 20 $file > ${file/.sorted.bam}.filtered20.bam && \
samtools sort -O bam -o ${file/.sorted.bam}.filtered20.sorted.bam ${file/.sorted.bam}.filtered20.bam && \
samtools index -c ${file/.sorted.bam}.filtered20.sorted.bam" >> filter; done
ls5_launcher_creator.py -j filter -n filter -t 00:20:00 -N 1 -w 4 -a mega2014 -e cbscott@utexas.edu -q development

for file in *filtered20.sorted.bam; do
cat Coral_scaffolds | tr "\n" " " | xargs samtools view -bh $file > ${file/.sorted.bam}.coral.bam && echo $file && samtools view -c ${file/.sorted.bam}.coral.bam; done

for file in *filtered20.sorted.bam; do
samtools view -bh $file chr13 > ${file/.sorted.bam}.cladec.bam && echo $file && samtools view -c ${file/.sorted.bam}.cladec.bam; done

for file in *filtered20.sorted.bam; do
samtools view -bh $file chr11 > ${file/.sorted.bam}.cladea.bam && echo $file && samtools view -c ${file/.sorted.bam}.cladea.bam; done

#FILTER by second qual: 30
for file in *default.sorted.bam; do
echo " samtools view -bSq 30 $file > ${file/.sorted.bam}.filtered30.bam && \
samtools sort -O bam -o ${file/.sorted.bam}.filtered30.sorted.bam ${file/.sorted.bam}.filtered30.bam && \
samtools index -c ${file/.sorted.bam}.filtered30.sorted.bam" >> filter; done
ls5_launcher_creator.py -j filter -n filter -t 00:05:00 -N 1 -w 4 -a tagmap -e cbscott@utexas.edu -q development

for file in *filtered30.sorted.bam; do
cat Coral_scaffolds | tr "\n" " " | xargs samtools view -bh $file > ${file/.sorted.bam}.coral.bam && echo $file && samtools view -c ${file/.sorted.bam}.coral.bam; done

for file in *filtered30.sorted.bam; do
samtools view -bh $file chr13 > ${file/.sorted.bam}.cladec.bam && echo $file && samtools view -c ${file/.sorted.bam}.cladec.bam; done

for file in *filtered30.sorted.bam; do
samtools view -bh $file chr11 > ${file/.sorted.bam}.cladea.bam && echo $file && samtools view -c ${file/.sorted.bam}.cladea.bam; done

#index and sort the segregated bams
ls *.coral.bam >> subsets
ls *.cladea.bam >> subsets
ls *.cladec.bam >> subsets

for file in `cat subsets`; do
samtools sort -O bam -o ${file/.bam}.sorted.bam ${file/.bam}.bam && \
samtools index -c ${file/.bam}.sorted.bam; done

#MAPDAM for all, manually evaulate profiles
module load gsl
>total_mapdam
for file in *.sorted.bam; do
echo "mapDamage --input $file -r $GENOME_FASTA --merge-libraries --no-stats" >> total_mapdam; done
ls5_launcher_creator.py -j total_mapdam -n total_mapdam -t 00:15:00 -N 1 -w 4 -a tagmap -e cbscott@utexas.edu -q development

for file in *.sorted.bam; do
mapDamage --input $file -r $GENOME_FASTA --merge-libraries --no-stats; done


#What if the other reads didn't show promise due to the mapping method? What if there was ancient coral in the onshore cores?
>convert
for file in *.bam; do
echo "samtools fastq $file > ${file/.bam}.fastq" >> convert; done
ls5_launcher_creator.py -j convert -n convert -t 00:05:00 -w 4 -N 1 -a tagmap -e cbscott@utexas.edu -q development
sbatch convert.slurm
rm convert* #keep this directory cleaner

export GENOME_FASTA=$WORK/ancientDNAdb/Apalm_and_sym.fasta #have indexed ref genome here. Be careful, it seems we accidentally mapped to the WRONG genome off the bat.
>maps_default
for file in *.fastq; do
echo "bowtie2 --no-unal -x $GENOME_FASTA -U $file -S ${file/.fastq}.default.sam && \
samtools sort -O bam -o ${file/.fastq/}.default.sorted.bam ${file/.fastq/}.default.sam && samtools index -c ${file/.fastq/}.default.sorted.bam " >> maps_default;
done
ls5_launcher_creator.py -j maps_default -n maps_default -t 00:45:00 -w 4 -N 1 -a tagmap -e cbscott@utexas.edu -q normal
sbatch maps_default.slurm #these reads do not show promise


####EXPlORE EFFICIENCY OF PMD tools
export file=S17463.default.sorted.bam
echo "samtools view -h $file | pmdtools.0.60.py --threshold 0.5 --header --UDGhalf | samtools view -Sb - > ${file/.bam}.pmds0_5filter.bam" > pmd_filter
ls5_launcher_creator.py -j pmd_filter -n pmd_filter -t 00:15:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q development

export file=S17463.default.cladec.sorted.bam
samtools view -h $file | pmdtools.0.60.py -p --UDGhalf > udg_pmdscores.txt

for file in *.bam; do
echo "samtools sort -O bam -o ${file/.bam}.sorted.bam ${file/.bam}.bam && \
samtools index -c ${file/.bam}.sorted.bam" >> filter; done
ls5_launcher_creator.py -j filter -n filter -t 00:10:00 -N 1 -w 2 -a tagmap -e cbscott@utexas.edu -q development

for file in *.sorted.bam; do
cat Coral_scaffolds | tr "\n" " " | xargs samtools view -bh $file > ${file/.sorted.bam}.coral.bam && echo $file && samtools view -c ${file/.sorted.bam}.coral.bam; done

for file in *.sorted.bam; do
samtools view -bh $file chr13 > ${file/.sorted.bam}.cladec.bam && echo $file && samtools view -c ${file/.sorted.bam}.cladec.bam; done

for file in *.sorted.bam; do
samtools view -bh $file chr11 > ${file/.sorted.bam}.cladea.bam && echo $file && samtools view -c ${file/.sorted.bam}.cladea.bam; done

samtools view -bSq 30 S17463.default.sorted.pmds1_5filter.cladec.bam | samtools view -c


###CLIPPED MAPPING ANALYSIS- NEED TO REDO S17465
#clipping is still a major issue
for file in *clip.sorted.bam; do
cat ../map2_coral_zoox/Coral_scaffolds | tr "\n" " " | xargs samtools view -bh $file > ${file/.sorted.bam}.coral.bam && echo $file && samtools view -c ${file/.sorted.bam}.coral.bam; done

for file in *clip.coral.bam; do
  samtools sort -O bam -o ${file/.bam}.sorted.bam ${file/.bam}.bam && \
  samtools index -c ${file/.bam}.sorted.bam; done

for file in *.sorted.bam; do
mapDamage --input $file -r $GENOME_FASTA --merge-libraries --no-stats; done

### Deal with stupid naming of Voolstra Genomes
awk -F',' '{print $2}' id_lookup | awk -F'/' '{print $6}' > nonsense_ids_ord

for f in nonsense_ids_ord sense_ids_ord; do sed -i "s/$/\t$f/" $f; done

#Gotta trim fastas for known damage, then remap.
#create new fastas from the coral reads used for angsd ONLY
sed 's/^[ATCG]\(.*\).$/\1/' S17463.fastq > test.fastq
pwd = /scratch/06909/cbscott/ancientDNA/aDNA_check/compareWGS/apalm_map

for file in S174*.bam; do
samtools fastq $file > ${file/.sorted.bam}.fastq; done
# now clip first and last base
module load cutadapt
for file in S174*.coral.fastq; do
cutadapt -u 1 -o ${file/.fastq}.temp.fastq $file; done
for file in S174*.temp.fastq; do
cutadapt -u -1 -o ${file/.temp.fastq}.clip.fastq $file && rm $file; done

#remap
>maps_default
for file in *.clip.fastq; do
echo "bowtie2 --no-unal -x $GENOME_FASTA -U $file -S ${file/.ancient.coral.clip.fastq/}.clip.coral.sam && \
samtools sort -O bam -o ${file/.ancient.coral.clip.fastq/}.clip.coral.sorted.bam ${file/.ancient.coral.clip.fastq/}.clip.coral.sam && samtools index -c ${file/.ancient.coral.clip.fastq/}.clip.coral.sorted.bam  " >> maps_default;
done
ls5_launcher_creator.py -j maps_default -n maps_default -t 00:10:00 -N 1 -w 2 -a mega2014 -e cbscott@utexas.edu -q development

#Refilter for coral reads
for file in *.clip.coral.sorted.bam; do
samtools view -bh $file `cat Coral_scaffolds` > ${file/.coral.sorted.bam}.coral.bam && echo $file && samtools view -c ${file/.coral.sorted.bam}.coral.bam; done

for file in *.clip.coral.bam; do
  samtools sort -O bam -o ${file/.bam}.sorted.bam ${file/.bam}.bam && \
  samtools index -c ${file/.bam}.sorted.bam; done
#done with clipped smalls

#IS it cladicopium? Get new NCBI ref.

$HOME/datasets download assembly GCA_003297045.1 #SymC
$HOME/datasets download assembly GCA_003297005.1 #SymA
$HOME/datasets download assembly GCA_000507305.1 #SymB
#no ncbi dataset for D, so
sed -n '/>chr14/,$p' $WORK/ancientDNAdb/symABCD_genome.fasta > symD_ref.fna


export SYMREF=$WORK/symref/ncbi_syma_symc_genomic.fna

echo "bowtie2-build $SYMREF $SYMREF" >btb
ls5_launcher_creator.py -j btb -n btb -l btbl -t 01:00:00 -a mega2014 -e cbscott@utexas.edu -w 1 -q normal
sbatch btbl
samtools faidx $SYMREF

>maps_default
for file in *.fastq; do
echo "bowtie2 --no-unal -x $SYMREF -U $file -S ${file/.fastq}.newsym.sam && \
samtools sort -O bam -o ${file/.fastq/}.newsym.sorted.bam ${file/.fastq/}.newsym.sam && samtools index -c ${file/.fastq/}.newsym.sorted.bam " >> maps_default;
done
ls5_launcher_creator.py -j maps_default -n maps_default -t 00:45:00 -w 4 -N 1 -a mega2014 -e cbscott@utexas.edu -q normal
sbatch maps_default.slurm


export A_REF=/work/06909/cbscott/lonestar/symref/GCA_003297005.1_SymA_ver_1.0_genomic.fna
export C_REF=/work/06909/cbscott/lonestar/symref/GCA_003297045.1_SymC_ver_1.0_genomic.fna

grep ">" $A_REF | sed 's/^.//' | awk '{print $1}' > a_scaffolds
grep ">" $C_REF | sed 's/^.//' | awk '{print $1}' > c_scaffolds

for file in *.sorted.bam; do
cat a_uniq_scaffolds | tr "\n" " " | xargs samtools view -bh $YOURBAM > $YOURBAM.subset && echo $file && samtools view -c ${file/.sorted.bam}.a.bam; done
#doesn't work.
for file in *.sorted.bam; do
samtools view -bh $file `cat a_scaffolds` > ${file/.sorted.bam}.a.bam && echo $file && samtools view -c ${file/.sorted.bam}.a.bam; done #works

for file in *.sorted.bam; do
samtools view -bh $file `cat c_scaffolds` > ${file/.sorted.bam}.c.bam && echo $file && samtools view -c ${file/.sorted.bam}.c.bam; done #works

#Now remapdamage for these reads.
export SYMREF=$WORK/symref/ncbi_syma_symc_genomic.fna

for file in *c.bam; do
  samtools sort -O bam -o ${file/.bam}.sorted.bam ${file/.bam}.bam && \
  samtools index -c ${file/.bam}.sorted.bam; done


  #MAPDAM for all, manually evaulate profiles
module load gsl
>total_mapdam
for file in *.sorted.bam; do
echo "mapDamage --input $file -r $SYMREF --merge-libraries --no-stats" >> total_mapdam; done
ls5_launcher_creator.py -j total_mapdam -n total_mapdam -t 00:15:00 -N 1 -w 1 -a mega2014 -e cbscott@utexas.edu -q development

#filter and remap just to check
for file in *.sorted.bam; do
samtools view -bSq 20 $file > ${file/.sorted.bam}.filtered20.bam && \
samtools sort -O bam -o ${file/.sorted.bam}.filtered20.sorted.bam ${file/.sorted.bam}.filtered20.bam && \
samtools index -c ${file/.sorted.bam}.filtered20.sorted.bam; done

for file in *.filtered20.sorted.bam; do
mapDamage --input $file -r $SYMREF --merge-libraries --no-stats; done
ls5_launcher_creator.py -j total_mapdam -n total_mapdam -t 00:15:00 -N 1 -w 1 -a mega2014 -e cbscott@utexas.edu -q development


#CREATE MASTER GENOME
export MASTER_GEN=$WORK/ancientDNAdb/mastergenome/master_mixed_genome.fasta
echo "bowtie2-build $MASTER_GEN $MASTER_GEN" >btb
ls5_launcher_creator.py -j btb -n btb -l btbl -t 04:00:00 -a mega2014 -e cbscott@utexas.edu -w 1 -q normal
sbatch btbl
samtools faidx $MASTER_GEN


>maps_default
for file in *.fastq; do
echo "bowtie2 --no-unal -x $MASTER_GEN -U $file -S ${file/.fastq}.master.sam && \
samtools sort -O bam -o ${file/.fastq/}.master.sorted.bam ${file/.fastq/}.master.sam && samtools index -c ${file/.fastq/}.master.sorted.bam " >> maps_default;
done
ls5_launcher_creator.py -j maps_default -n maps_default -t 04:00:00 -a mega2014 -e cbscott@utexas.edu -w 4 -q normal

#Now need to separate all of these reads by organism...
export A_REF=/work/06909/cbscott/lonestar/symref/GCA_003297005.1_SymA_ver_1.0_genomic.fna
export C_REF=/work/06909/cbscott/lonestar/symref/GCA_003297045.1_SymC_ver_1.0_genomic.fna
export B_REF=$WORK/symref/ref_gen_B.fna
grep ">" $A_REF | sed 's/^.//' | awk '{print $1}' > a_scaffolds
grep ">" $B_REF | sed 's/^.//' | awk '{print $1}' > b_scaffolds
grep ">" $C_REF | sed 's/^.//' | awk '{print $1}' > c_scaffolds

#Filter by qual first
for file in *.sorted.bam; do
samtools view -bSq 20 $file > ${file/.sorted.bam}.filtered20.bam && \
samtools sort -O bam -o ${file/.sorted.bam}.filtered20.sorted.bam ${file/.sorted.bam}.filtered20.bam && \
samtools index -c ${file/.sorted.bam}.filtered20.sorted.bam; done
ls5_launcher_creator.py -j filter -n filter -t 00:20:00 -N 1 -w 4 -a mega2014 -e cbscott@utexas.edu -q development

>names
>counts
for file in *.master.filtered20.sorted.bam; do
samtools view -bh $file `cat Coral_scaffolds` > ${file/.master.filtered20.sorted.bam}.coral.filtered20.bam && echo ${file/.master.filtered20.sorted.bam}.coral.filtered20.bam >> names && samtools view -c ${file/.master.filtered20.sorted.bam}.coral.filtered20.bam >> counts && \
samtools view -bh $file `cat a_scaffolds` > ${file/.master.filtered20.sorted.bam}.a.filtered20.bam && echo ${file/.master.filtered20.sorted.bam}.a.filtered20.bam >> names && samtools view -c ${file/.master.filtered20.sorted.bam}.a.filtered20.bam >> counts && \
samtools view -bh $file `cat b_scaffolds` > ${file/.master.filtered20.sorted.bam}.b.filtered20.bam && echo ${file/.master.filtered20.sorted.bam}.b.filtered20.bam >> names && samtools view -c ${file/.master.filtered20.sorted.bam}.b.filtered20.bam >> counts && \
samtools view -bh $file `cat c_scaffolds` > ${file/.master.filtered20.sorted.bam}.c.filtered20.bam && echo ${file/.master.filtered20.sorted.bam}.c.filtered20.bam >> names && samtools view -c ${file/.master.filtered20.sorted.bam}.c.filtered20.bam >> counts && \
samtools view -bh $file chr14 > ${file/.master.filtered20.sorted.bam}.d.filtered20.bam && echo ${file/.master.filtered20.sorted.bam}.d.filtered20.bam >> names && samtools view -c ${file/.master.filtered20.sorted.bam}.d.filtered20.bam >> counts && \
samtools view -bh $file `cat cca_scaffolds` > ${file/.master.filtered20.sorted.bam}.cca.filtered20.bam && echo ${file/.master.filtered20.sorted.bam}.cca.filtered20.bam >> names && samtools view -c ${file/.master.filtered20.sorted.bam}.cca.filtered20.bam >> counts && \
samtools view -bh $file `cat red_scaffolds` > ${file/.master.filtered20.sorted.bam}.red.filtered20.bam && echo ${file/.master.filtered20.sorted.bam}.red.filtered20.bam >> names && samtools view -c ${file/.master.filtered20.sorted.bam}.red.filtered20.bam >> counts; done

#need to create a file that doesn't have any of these scaffolds OH can just pull them from the og metagenome....
export METAREF=$SCRATCH/ancientDNA/Voolstra_Genomes/proper_name/labeled_metaref.fasta
grep ">" $METAREF | sed 's/^.//' | awk '{print $1}' > meta_scaffolds

for file in *.master.filtered20.sorted.bam; do
samtools view -bh $file `cat meta_scaffolds` > ${file/.master.filtered20.sorted.bam}.meta.filtered20.bam && echo ${file/.master.filtered20.sorted.bam}.meta.filtered20.bam >> names && samtools view -c ${file/.master.filtered20.sorted.bam}.meta.filtered20.bam >> counts; done

#all the filtered 20s now need to be sorted and indexed....
>sort
for file in *.filtered20.bam; do
echo "samtools sort -O bam -o ${file/.filtered20.bam}.filtered20.sorted.bam $file && \
samtools index -c ${file/.filtered20.bam}.filtered20.sorted.bam" >> sort; done
ls5_launcher_creator.py -j sort -n sort -t 00:20:00 -N 1 -w 1 -a mega2014 -e cbscott@utexas.edu -q development


export MASTER_GEN=$WORK/ancientDNAdb/mastergenome/master_mixed_genome.fasta
metAchecker_1.sh S17464.meta.filtered20.sorted.bam $MASTER_GEN "S17464"




#MODIFIED aCHECKER for coral, zoox, reds, cca....
#create list of relevant files
>nonmeta_check
ls *.coral*sorted.bam >> nonmeta_check
ls *.c.*sorted.bam >> nonmeta_check
ls *.a.*sorted.bam >> nonmeta_check
ls *.b.*sorted.bam >> nonmeta_check
ls *.d.*sorted.bam >> nonmeta_check
ls *.cca.*sorted.bam >> nonmeta_check
ls *.red.*sorted.bam >> nonmeta_check
ls *.meta.filtered20.sorted.bam >> nonmeta_check

for file in `cat nonmeta_check`; do
mapDamage --input $file -r $WORK/ancientDNAdb/mastergenome/master_mixed_genome.fasta --merge-libraries --no-stats; done

for file in *.sorted.mapDamage; do
metAchecker_likelyancient.R $file ${file/.sorted.mapDamage}.tempo; done #might still need to be a job?

for file in *.tempo; do
awk -F' ' '{print $2}' $file > ${file/.tempo}.check && rm $file; done

mkdir overall.ancient_mapping
mkdir overall.modern_mapping

>overall.likely_ancient
>overall.likely_modern
for file in *.check; do
metA_remover.sh $file overall && rm $file; done

S17464.meta.filtered20

metAchecker_tables.R S17464.meta.filtered20.aligns S17464.meta.filtered20.full_names S17464.meta.filtered20.top_hits
awk -F' ' '{print $2}' S17464.meta.filtered20.top_hits | sed 's/\"//g' > S17464.meta.filtered20.alt_tops

#PIA_metagenomics???
#In theory, don't have to build a db, but can use already downloaded blast nt db
#Still have to build it?
export NT_DB=$WORK/update_ntdb/nt
export TEST_FA=/scratch/06909/cbscott/ancientDNA/aDNA_check/PIA_metagenome/S17463_test_named.fa

#make blast db
echo "makeblastdb -in nt.fa -dbtype nucl" > largemem_db
module load TACC-largemem
ls5_launcher_creator.py -j largemem_db -n largemem_db -t 04:00:00 -N 1 -w 1 -a mega2014 -e cbscott@utexas.edu -q largemem512GB
#if this doesn't work, also creating a database direct from NCBI's update_db command
#NEED fastqs as fasta.. should be job... these are large files.
echo "paste - - - - < $TEST_FA | cut -f 2 | sed 's/^@/>/' | tr '\t' '\n' > S17463_noname.fa" > fa_convert
echo "paste - - - - < $TEST_FA | cut -f 2 | sed 's/^@/>/' | tr '\n' '\t' > S17463_noname.fa" > fa_convert

blastn –db $NT_DB -num_threads 1 –query $TEST_FA –out ${TEST_FA/.fa}.out  -max_target_seqs 500 -outfmt "6 std staxids"

module load TACC-largemem
ls5_launcher_creator.py -j largmem_short -n largmem_short -t 00:10:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q development

##Sequence length versus percent damaged plot...
#Break bams into size classes:
export file=S17463.coral.filtered20.sorted.bam
for file in *.sorted.bam; do
samtools view -h $file | awk 'length($10) <= 29 || $1 ~ /^@/' | samtools view -bS - > ${file/.sorted.bam}.30minus.size.bam && \
samtools view -h $file | awk 'length($10) > 29 && length($10) < 40 || $1 ~ /^@/' | samtools view -bS - > ${file/.sorted.bam}.30t40.size.bam && \
samtools view -h $file | awk 'length($10) > 39 && length($10) < 50 || $1 ~ /^@/' | samtools view -bS - > ${file/.sorted.bam}.40t50.size.bam && \
samtools view -h $file | awk 'length($10) > 49 && length($10) < 60 || $1 ~ /^@/' | samtools view -bS - > ${file/.sorted.bam}.50t60.size.bam && \
samtools view -h $file | awk 'length($10) > 59 && length($10) < 70 || $1 ~ /^@/' | samtools view -bS - > ${file/.sorted.bam}.60t70.size.bam && \
samtools view -h $file | awk 'length($10) > 69 && length($10) < 80 || $1 ~ /^@/' | samtools view -bS - > ${file/.sorted.bam}.70t80.size.bam && \
samtools view -h $file | awk 'length($10) > 79 && length($10) < 90 || $1 ~ /^@/' | samtools view -bS - > ${file/.sorted.bam}.80t90.size.bam && \
samtools view -h $file | awk 'length($10) >= 90 || $1 ~ /^@/' | samtools view -bS - > ${file/.sorted.bam}.90plus.size.bam; done

#do a little mapDam test first


export MASTER_GEN=$WORK/ancientDNAdb/mastergenome/master_mixed_genome.fasta
module load gsl

>total_mapdam
for file in *.size.bam; do
echo "mapDamage --input $file -r $MASTER_GEN --merge-libraries --no-stats" >> total_mapdam; done
ls5_launcher_creator.py -j total_mapdam -n total_mapdam -t 00:15:00 -N 1 -w 4 -a tagmap -e cbscott@utexas.edu -q development
