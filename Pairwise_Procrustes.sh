##Pairwise Procrustes

######======CREATE SITES/SCAFFOLDS LISTS=======########
#start in location with all bams, then move things into a more specific directory
export file="S17466" #this is your ancient filename
echo ${file}.softclip.sorted.coral.bam > S17466
JUST_ANC="-minMapQ 20 -minQ 20 -doMaf 1 -doMajorMinor 1 -GL 1 -out anc_${file}"
angsd -b $file $JUST_ANC


gunzip anc_${file}.mafs.gz
cut -f 1,2,3,4 anc_${file}.mafs | tail -n +2 > sites_${file}.txt #need only first four columns, without a header
#cut -f1 sites_${file}.txt | sort | uniq >chrs_${file}.txt
angsd sites index sites_${file}.txt

#### move sites_*.txt and chrs_*.txt locally for Rscript:
#Run scaffold_sampler.R

for f in sampled_sites_${file}*.txt; do
  cut -d' ' -f1,4,5,6 $f > formatted_${f}; done
#this needs to be a job in the future.
>index
for f in formatted_sampled_sites_${file}*.txt; do
  echo "cut -d' ' -f1 $f | uniq > formatted_sampled_chrs_${f#formatted_sampled_sites_} && angsd sites index $f" >> index; done
  ls5_launcher_creator.py -j index -n index -a tagmap -e cbscott@utexas.edu -t 00:05:00 -w 10 -N 1 -q development
sbatch index.slurm
#Move output from scaffold_sampler back to TACC

>chrmaker
for file in sampled_sites_S17463*.txt; do
  echo "awk '{print \$1}' $file | sort | uniq > sampled_chrs_${file#sampled_sites_}" >> chrmaker; done #no longer need this, doing it in the indexing step
  ls5_launcher_creator.py -j chrmaker -n chrmaker -a tagmap -e cbscott@utexas.edu -t 00:30:00 -w 10 -N 1 -q development
#####=======RUN ANGSD REPLICATES==========##########
export MinIndPerc=0.5
export file="S17466"
export bams="bams_17466"
export out="myresult17466"

FILTERS0="-minInd $MI -minMapQ 20 -minQ 20 -minMaf 0.039 -dosnpstat 1 -doHWE 1 -hetbias_pval 1e-3 -skipTriallelic 1"
echo 'export NIND=`cat $bams | wc -l`; export MI=`echo "($NIND*$MinIndPerc+0.5)/1" | bc`' >calc1
for i in {1..100}; do
echo "source calc1 && angsd -b $bams -GL 1 $FILTERS0 -doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doPost 1 -doGlf 2 -doQsDist 1 -sites formatted_sampled_sites_${file}_${i}.txt -rf formatted_sampled_chrs_${file}_${i}.txt -P 12 -out ${out}.${i} && Rscript ~/bin/detect_clones.R $bams ${out}.${i}.ibsMat 0.15">>a_${file};
done
ls5_launcher_creator.py -j a_${file} -n a_${file} -a tagmap -e cbscott@utexas.edu -t 03:00:00 -w 16 -N 4 -q normal
sbatch a_${file}.slurm


#Move IBS matrices locally & run resampleIBS.R
