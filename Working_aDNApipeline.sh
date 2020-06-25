#AncientDNA things that actually work :)
#####################################################
#########     GET READS       #######################
#####################################################

#####################################################
#########    BUILD DATABASE   #######################
#####################################################
#plantless-malt contains invertebrate, bacteria, archaea, protista, and ALGAL genomes (even though algae are technically a plant- trying to avoid things like oak trees here)
echo "malt-build --input $SCRATCH/ancientDNA/bact_refgenomes/plantless-malt/* -s DNA --index $SCRATCH/ancientDNA/malt-allmarine --acc2taxonomy $SCRATCH/ancientDNA/bact_refgenomes/NCBI_kraken/taxonomy/nucl_all.accession2taxid --threads 24 --verbose" > malt_noplant
ls5_launcher_creator.py -j malt_noplant -n malt_noplant -t 02:00:00 -N 1 -w 24 -a tagmap -e cbscott@utexas.edu -q development





#####################################################
#########     MAP READS TO DB    #######################
#####################################################
##### 1. Divide large fastq files into smaller files to avoid Java memory errors
split -l 500000 S17464.Y1.E1.L1.fastq 464small  #probably only need 20 parts, not 450, based on memory usage
ls 464smalla* > 464smallfiles

##### 2. Run malt on each subset of the file.
> 464_all
for file in `cat 464smallfiles`;
do echo "malt-run -i $SCRATCH/ancientDNA/rawreads/${file} -a $SCRATCH/ancientDNA/malt_tests --index $SCRATCH/ancientDNA/malt-allmarine/ --mode BlastN --format Text --alignmentType SemiGlobal" >> 464_all;
done
ls5_launcher_creator.py -j 464_all -n 464_all -t 04:00:00 -a tagmap -e cbscott@utexas.edu -N 10 -w 1 -q normal

#still encountering memroy errors. Try a very tiny file.
###### 3. Combine the outputs back together
cat 464small*.blastn.gz > 464_all.blastn.gz #this will take a couple minutes

###### 4. Visualize the output. How many reads mapped? What was in the sample? Make sure the settings are 'loose', we don't expect much to map to interesting species #use sedaDNA params
#ms=50
#me=0.01
#supp=0
#mpi=95 from Ambrechet et al

#will 10 minutes be enough time? It is for small files, for the large file, needed much longer. Took 20 minutes for full file
echo "blast2rma -i 464_assembled.blastn.gz -r $SCRATCH/ancientDNA/rawreads/S17464.Y1.E1.L1.fastq -f BlastText -o /scratch/06909/cbscott/ancientDNA/malt_tests -supp 0 -ms 50 -me 0.01 -mpi 95" > 464_assembledVIZ
ls5_launcher_creator.py -j 464_assembledVIZ -n 464_assembledVIZ -t 01:30:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q development

#---------I think the database is missing algae... this would improve the mapping of the unidentified reads?
#####################################################
######### FILTER READS FOR ANCIENT DNA  ##############
#####################################################
#Use PMDTools here


#####################################################
#########     MAP READS TO REFERENCE   ##############
#####################################################
/home1/06909/cbscott/bam2fastq-1.1.0/bam2fastq S17463.Y1.E1.L1.pmds3filter.bam -o S17463.ancientDNA.fastq #make a fastq file input

echo "bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x /work/06909/cbscott/lonestar/ancientDNAdb/Apalm_sym_euk.fasta -U S17463.ancientDNA.fastq -S S17463.Y1.E1.L1.sam && samtools sort -O bam -o S17463.Y1.E1.L1.sorted.bam S17463.Y1.E1.L1.sam && samtools index -c S17463.Y1.E1.L1.sorted.bam" > map_aReads
ls5_launcher_creator.py -j map_aReads -n map_aReads -t 00:20:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q development
