#Need to coerce the genome to look like the human genome....

#Create 24 chromosomes & remap (just to this, not zoox)
concatFasta.pl fasta=Apalm_assembly_v2.0_180910.scaffolds.fasta num=24
export GEN_RED=$WORK/ancientDNAdb/mastergenome/human-fied_genome/Apalm_assembly_v2_cc.fasta

#chromosomes need to be just numbers....
awk '/^>/{print ">" ++i; next}{print}' < $GEN_RED > numbered_reduced_Apalmata.fasta
export GEN_NEW=$WORK/ancientDNAdb/mastergenome/human-fied_genome/numbered_reduced_Apalmata.fasta

#build indexed genome from genred
echo "bowtie2-build $GEN_NEW $GEN_NEW" >btb
ls5_launcher_creator.py -j btb -n btb -l btbl -t 01:00:00 -a tagmap -e cbscott@utexas.edu -w 1 -q normal
sbatch btbl

samtools faidx $GEN_NEW
