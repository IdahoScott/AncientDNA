#Get aDNA files
wget https://www.dropbox.com/sh/vl2vwz1yc9kyc9u/AAAVRQ9RJgAJhzyXQ2IUuxtXa?dl=0 #rename file as a .zip, then extract

#add to .bin a path to bam2fastq
>tofastq
for file in *.bam; do
echo "/home1/06909/cbscott/bam2fastq-1.1.0/bam2fastq $file -o ${file/.bam/}.fastq" >> tofastq;
done

#Assumed reads were already trimmed
>maps
for file in *.fastq; do
echo "bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x $GENOME_FASTA -U $file -S ${file/.fastq/}.sam && \
samtools sort -O bam -o ${file/.fastq/}.sorted.bam ${file/.fastq/}.sam && samtools index ${file/.fastq/}.bam " >> maps;
done

ls5_launcher_creator.py -j maps_c -n maps_c -t 2:00:00 -w 24 -a tagmap -e cbscott@utexas.edu -q development


#For .csi indexes. This will be necessary because the massive size of the reference genome
>maps_c
for file in S17463.Y1.E1.L1.fastq; do
echo "bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x $GENOME_FASTA -U $file -S ${file/.fastq/}.sam && \
samtools sort -O bam -o ${file/.fastq/}.sorted.bam ${file/.fastq/}.sam && samtools index -c ${file/.fastq/}.sorted.bam " >> maps_c;
done

for file in *.bam; do
echo "samtools index -c $file" >> maps_c;
done

#---------------------------From Ambrechet et al 2020
#MEGAN6 (version 6_15_1; http://ab.inf.uni-tuebingen.de/software/megan6/; required to run Blast2RMA)
#wget https://software-ab.informatik.uni-tuebingen.de/download/megan6/MEGAN_Community_unix_6_18_10.sh
#chmod +x MEGAN_Community_unix_6_18_10.sh
#sh MEGAN_Community_unix_6_18_10.sh
#add path in .bashrc

#Failed?
#wget https://software-ab.informatik.uni-tuebingen.de/download/malt/MALT_unix_0_4_1.sh
#chmod +x MALT_unix_0_4_1.sh
#sh MALT_unix_0_4_1.sh
#Actually, probably don't need new programs from Ambrechet - just create a bigger genome with SILVA database
#Why SILVA?
#SILVA SSU database was chosen due to the extensive use of the 18S gene in the marine eukaryote taxonomy literature

#PARALLEL
cd
wget http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2
tar -xvf parallel-latest.tar.bz2
#From the README
(wget -O - pi.dk/3 || lynx -source pi.dk/3 || curl pi.dk/3/ || \
     fetch -o - http://pi.dk/3 ) > install.sh
sha1sum install.sh | grep 3374ec53bacb199b245af2dda86df6c912345678374ec53bacb199b245af2dda86df6c9
md5sum install.sh | grep 029a9ac06e8b5bc6052eac57b2c3c9ca029a9ac06e8b5bc6052eac57b2c3c9ca
sha512sum install.sh | grep f517006d9897747bed8a4694b1acba1b40f53af69e20dae5713ba06cf517006d9897747bed8a4694b1acba1b1464beb4600556293f2356f33e9c4e3c76e3f3afa9db4b32bd33322b975696fce6b23cfb
bash install.sh


#try using MALT? #Took ~12 minutes
echo "malt-build -i $WORK/ancientDNAdb/SILVA_132_SSURef_Nr99_tax_silva.fasta -s DNA -d /work/06909/cbscott/lonestar/ancientDNAdb/malt-index/ --threads 16 --verbose" > buildmalt
ls5_launcher_creator.py -j buildmalt -n buildmalt -t 04:00:00 -a tagmap -e cbscott@utexas.edu -N 1 -w 16 -q normal #only took ten minutes?
#this works, now run MALT, just try for one file
echo "malt-run -i $SCRATCH/ancientDNA/rawreads/S17463.Y1.E1.L1.fastq -a $SCRATCH/ancientDNA/rawreads/ --index $WORK/ancientDNAdb/malt-index/ --mode BlastN --format Text --alignmentType SemiGlobal -t 16 --verbose"> runmalt
#Return a .sam file
echo "malt-run -i $SCRATCH/ancientDNA/rawreads/test.fastq -a $SCRATCH/ancientDNA/malt_tests --index $SCRATCH/ancientDNA/malt-allmarine/ --mode BlastN --format Tab --alignmentType SemiGlobal" > malt_sub
ls5_launcher_creator.py -j malt_sub -n malt_sub -t 00:15:00 -a tagmap -e cbscott@utexas.edu -N 1 -w 1 -q development

echo "malt-run -i $SCRATCH/ancientDNA/rawreads/test.fastq -a $SCRATCH/ancientDNA/malt_tests --index $SCRATCH/ancientDNA/malt-allmarine/ --mode BlastN --format Text --alignmentType SemiGlobal" > malt_sub
ls5_launcher_creator.py -j malt_sub -n malt_sub -t 00:05:00 -a tagmap -e cbscott@utexas.edu -N 1 -w 1 -q development

#Try with whole file again? no threading
echo "malt-run -i $SCRATCH/ancientDNA/rawreads/464.fastq -a $SCRATCH/ancientDNA/malt_tests --index $SCRATCH/ancientDNA/malt-allmarine/ --mode BlastN --format Text --alignmentType SemiGlobal" > malt_all
ls5_launcher_creator.py -j malt_all -n malt_all -t 01:30:00 -a tagmap -e cbscott@utexas.edu -N 1 -w 1 -q normal #is development causing the issue here?


#Try to automate entire file
> 464_all
for file in `cat 464files`;
do echo "malt-run -i $SCRATCH/ancientDNA/rawreads/${file} -a $SCRATCH/ancientDNA/malt_tests --index $SCRATCH/ancientDNA/malt-allmarine/ --mode BlastN --format Text --alignmentType SemiGlobal" >> 464_all;
done
ls5_launcher_creator.py -j 464_all -n 464_all -t 00:10:00 -a tagmap -e cbscott@utexas.edu -N 5 -w 1 -q normal

#Still too big.... reduce file size by 4

#Now run blast2RMA
#blast2RMA (requires lots of memory)
#Script for blast2RMA named run-blast2rma.sh
input=$SCRATCH/ancientDNA/rawreads/ #the directory with the de-dupe.gz files, in our case, the fastq file S17463.Y1.E1.L1.fastq
ms=50
me=0.01
supp=0
mpi=95
i=S17463.Y1.E1.L1.fastq
b=${i/.fastq/.blastn.gz} #CANNOT SEEM TO FIX JAVA ERRORS? EVEN AFTER REINSTALLING JAVA, STILL RETURNING ERRORS
echo "blast2rma -r $i -i $b -o $input --mode BlastN -v -ms $ms -me $me -f BlastText -supp $supp -mpi $mpi -alg weighted -lcp 80" > blast2rma
ls5_launcher_creator.py -j blast2rma -n blast2rma -t 02:00:00 -a tagmap -e cbscott@utexas.edu -N 1 -w 1 -q normal #only took ten minutes?
#Ran into java erros, attempted to resolve by reinstalling Java Runtime Environment, current version =
#From Ambrechet
#for i in $input/*.gz-dedupe.gz; do b=${i/.gz-dedupe.gz/.blastn.gz}; blast2rma -r $i -i $b -o $input -v -ms $ms -me $me -f BlastText -supp $supp -mpi $mpi -alg weighted -lcp 80; done | parallel -j 8 --delay 6
#run blast2RMA
nohup ./run-blast2rma.sh &

#try sam2rma?
i=S17463.Y1.E1.L1.fastq
b=${i/.fastq/.blastn.sam.gz} #CANNOT SEEM TO FIX JAVA ERRORS? EVEN AFTER REINSTALLING JAVA, STILL RETURNING ERRORS
echo "sam2rma -r $i -i $b -o $input --mode BlastN --propertiesFile Megan.def" > sam2rma
ls5_launcher_creator.py -j sam2rma -n sam2rma -t 02:00:00 -a tagmap -e cbscott@utexas.edu -N 1 -w 2 -q normal #only took ten minutes?





#Try to get reads from SILVA database, use the same DB as Ambrechet in order to compare results
wget https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_SSURef_Nr99_tax_silva.fasta.gz
gunzip SILVA_132_SSURef_Nr99_tax_silva.fasta.gz

#Create a massive reference genome of bacterial data + coral + syms
cp SILVA_132_SSURef_Nr99_tax_silva.fasta.gz SILVA_marine_euk.fasta
cat Apalm_assembly.fasta symABCD.fasta SILVA_marine_euk.fasta > Apalm_sym_euk.fasta

#this is taking forever - maybe limit data base?

echo "bowtie2-build $GENOME_FASTA $GENOME_FASTA" >btb
ls5_launcher_creator.py -j btb -n btb -l btbl -t 3:30:00 -a tagmap -e cbscott@utexas.edu -w 1 -q normal
sbatch btbl

export GENOME_FASTA=$WORK/ancientDNAdb/Apalm_sym_euk.fasta
samtools faidx $GENOME_FASTA

#DON'T NEED TO TRIM READS, JUST MAP
>tofastq
for file in *.bam; do
echo "/home1/06909/cbscott/bam2fastq-1.1.0/bam2fastq $file -o ${file/.bam/}.fastq" >> tofastq;
done

#For .csi indexes. This will be necessary because the massive size of the reference genome
>maps_c
for file in *.fastq; do
echo "bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x $GENOME_FASTA -U $file -S ${file/.fastq/}.sam && \
samtools sort -O bam -o ${file/.fastq/}.coralzoox.sorted.bam ${file/.fastq/}.sam && samtools index -c ${file/.fastq/}.coralzoox.sorted.bam " >> maps_c;
done
ls5_launcher_creator.py -j maps_c -n maps_c -t 8:00:00 -N 4 -w 1 -a tagmap -e cbscott@utexas.edu -q normal

#Ran into bowtie segmentation fault... large size?

#Let's try mapping to Marinobacter hydrocabonoclasticus (for fun)
rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/284/615/GCF_000284615.1_ASM28461v1/GCF_000284615.1_ASM28461v1_genomic.fna.gz .
export MARINE_BACT=$SCRATCH/ancientDNA/bact_refgenomes/Marinobacteria_HydroCarbono_genome.fna
echo "bowtie2-build $MARINE_BACT $MARINE_BACT" >btb
ls5_launcher_creator.py -j btb -n btb -l btbl -t 1:00:00 -a tagmap -e cbscott@utexas.edu -w 1 -q development
samtools faidx $MARINE_BACT

#Let's get MarineRef full genomes (~1000 individuals) #(may have accidentally grabbed protiens instead of organisms)
#This doesn't seem to work? Try a different way!
echo "wget https://s1.sfb.uit.no/public/mar/MarRef/BLAST/nucleotides/marref_nucleotides_V5.fna" > getMarineRef
ls5_launcher_creator.py -j getMarineRef -n getMarineRef -t 1:30:00 -a tagmap -e cbscott@utexas.edu -w 1 -q normal

echo "wget https://s1.sfb.uit.no/public/mar/MarRef/BLAST/assembly/marref_assembly_V5.fa" > MarineRefInds
ls5_launcher_creator.py -j MarineRefInds -n MarineRefInds -t 12:00:00 -a tagmap -e cbscott@utexas.edu -w 1 -q normal

# try a different method to get the MarineRef files you want. This way they can also be parallelized? Only 944 out of 970
#15 mins almost long enough, 1:30:00 should be mroe than enough
wget https://s1.sfb.uit.no/public/mar/Resources/kraken/dl_list_marref.txt
>getMarRef
for file in `cat dl_list_marref.txt`; do echo "wget $file" >> getMarRef ; done
ls5_launcher_creator.py -j getMarRef -n getMarRef -t 01:30:00 -a tagmap -e cbscott@utexas.edu -N 10 -w 24 -q normal


#just map one sample
echo "bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x $MARINE_BACT -U $file -S ${file/.fastq/}.sam && \
samtools sort -O bam -o ${file/.fastq/}.sorted.bam ${file/.fastq/}.sam && samtools index -c ${file/.fastq/}.sorted.bam " > maps_bact
ls5_launcher_creator.py -j maps_bact -n maps_bact -t 0:30:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q development
sbatch maps_bact

#Let's try using KRAKEN and BRAKEN to classify these reads!
#Download and build KRAKEN
cdh
wget http://github.com/DerrickWood/kraken2/archive/v2.0.8-beta.tar.gz
tar -vxf v2.0.8.beta.tar.gz
mkdir kraken2
cd v2.0.8-beta.tar.gz
./install_kraken2.sh $HOME/kraken2
#add relevant files to PATH, as prompted - NOT AS PROMPTED - do as usual

#Build default KRAKEN db (it will be like 100GB, according to website)
#bacterial, archaeal, and viral domains, along with the human genome and a collection of known vectors (UniVec_Core)
#We have multiple cores, so maybe we can paralellize this?
export DBNAME=$SCRATCH/ancientDNA/bact_refgenomes/CustomKrakenDB/
#May have skipped a crucial step
echo "kraken2-build --download-taxonomy --threads 12 --db $DBNAME --use-ftp" > gettax
ls5_launcher_creator.py -j gettax -n gettax -t 02:00:00 -N 1 -w 12 -a tagmap -e cbscott@utexas.edu -q normal
#can we thread it? can we use ftp to get it done slower?

echo "kraken2-build --standard --threads 24 --db $DBNAME" > buildKrakenDB #probably going to take a long time, try 24 threads?
ls5_launcher_creator.py -j buildKrakenDB -n buildKrakenDB -t 02:00:00 -N 2 -w 12 -a tagmap -e cbscott@utexas.edu -q normal
#Acutally try to make  a custom library using the reference sequences, by adding other things to database first
#Add sym genomes and coral genomes separately, especially after knowing taxid?
#Need to add a specific string to each chr
sed 's/\<chr11\>/&|kraken:taxid|256873/' symABCD_kraken.fasta > symABCD_kraken_edit.fasta #each sym has different taxid (A=256873, B=88618, C=34742, D=24806)
sed 's/\<chr12\>/&|kraken:taxid|88618/' symABCD_kraken_edit.fasta > symABCD_kraken_edit1.fasta #each sym has different taxid (A=256873, B=88618, C=34742, D=24806)
sed 's/\<chr13\>/&|kraken:taxid|35742/' symABCD_kraken_edit1.fasta > symABCD_kraken_edit2.fasta #each sym has different taxid (A=256873, B=88618, C=34742, D=24806)
sed 's/\<chr14\>/&|kraken:taxid|24806/' symABCD_kraken_edit2.fasta > symABCD_kraken_final.fasta #each sym has different taxid (A=256873, B=88618, C=34742, D=24806)

#attmpted to run on login node- might take too long- need to submit a job
echo "kraken2-build --add-to-library symABCD_kraken_final.fasta --db $DBNAME" > addsyms
ls5_launcher_creator.py -j addsyms -n addsyms -t 01:00:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q development


#Try to get all libraries we want
echo "kraken2-build --download-library 'archaea' --threads 24 --db $DBNAME --use-ftp" > getarchaea
ls5_launcher_creator.py -j getarchaea -n getarchaea -t 00:45:00 -N 2 -w 12 -a tagmap -e cbscott@utexas.edu -q normal
#finally! we can build our kraken DATABASE

echo "kraken2-build --db=$DBNAME --threads=24" > buildDB
ls5_launcher_creator.py -j buildDB -n buildDB -t 03:00:00 -N 2 -w 12 -a tagmap -e cbscott@utexas.edu -q normal #took 54 minutes
#no threads, should take ~1:30 let's get rid of the plants?

echo "kraken2-build --build --db=/scratch/06909/cbscott/ancientDNA/bact_refgenomes/NCBI_kraken/ --threads=1 --kmer-len 45" > buildsmallDB
ls5_launcher_creator.py -j buildsmallDB -n buildsmallDB -t 03:00:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q normal #took 1:30 minutes


#Now follow braken pipeline
export THREADS=24
export KRAKEN_DB=$SCRATCH/ancientDNA/bact_refgenomes/NCBI_kraken
echo "bracken-build -d ${KRAKEN_DB} -t ${THREADS}" > brakenbuild
ls5_launcher_creator.py -j brakenbuild -n brakenbuild -t 02:00:00 -N 2 -w 12 -a tagmap -e cbscott@utexas.edu -q normal


#kraken2-build --download-library bacteria --db $DBNAME
#kraken2-build --download-library human --db $DBNAME
#kraken2-build --download-library fungi --db $DBNAME
#kraken2-build --download-library plant --db $DBNAME
#kraken2-build --download-library protozoa --db $DBNAME" > get_libs
#ls5_launcher_creator.py -j get_libs -n get_libs -t 02:00:00 -N 2 -w 3 -a tagmap -e cbscott@utexas.edu -q normal

#add coral and syms to this library

echo "kraken2-build --add-to-library $WORK/ancientDNAdb/symABCD_kraken_edit.fasta --db $DBNAME" > addsyms
ls5_launcher_creator.py -j addsyms -n addsyms -t 01:00:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu


#Someone has already compiled a Kraken MAR db - this uses kraken1 not kraken2 -- that's why it's not working
wget https://s1.sfb.uit.no/public/mar/Resources/kraken/build-KrakenMar.sh
echo "./build-KrakenMar.sh" > KrakenMAR
ls5_launcher_creator.py -j KrakenMAR -n KrakenMAR -t 06:00:00 -N 2 -w 16 -a tagmap -e cbscott@utexas.edu -q normal

#jellyfish installation
#I think I may have been missing a dependency
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-2.3.0.tar.gz
tar -xvf jellyfish-2.3.0.tar.gz
cd jellyfish-2.3.0
./configure --PATH=$HOME
make
make install

#Libraries have been installed in:
#   /home1/06909/cbscott/lib

#If you ever happen to want to link against installed libraries
#in a given directory, LIBDIR, you must either use libtool, and
#specify the full pathname of the library, or use the '-LLIBDIR'
#flag during linking and do at least one of the following:
#   - add LIBDIR to the 'LD_LIBRARY_PATH' environment variable
#     during execution
  # - add LIBDIR to the 'LD_RUN_PATH' environment variable
#     during linking
#   - use the '-Wl,-rpath -Wl,LIBDIR' linker flag
#   - have your system administrator add LIBDIR to '/etc/ld.so.conf'
#PATH            -> /my/dir/bin
#MANPATH         -> /my/dir/share/man
#PKG_CONFIG_PATH -> /my/dir/lib/pkgconfig


#-----------GET KRAKEN DB BY HAND-----------
#Get NCBI taxonomy
export NEW_DB=$SCRATCH/ancientDNA/bact_refgenomes/NCBI_kraken
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz $NEW_DB/taxonomy #works
cd taxonomy
tar -xvf new_taxdump.tar.gz

mkdir library
cd library
mkdir archaea
echo "wget --recursive https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/ $NEW_DB/library/archaea" > getarch #doesn't work!!!!

#https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/ why won't ftp requests work?
#manually get all ftp links for a certain taxa
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt $NEW_DB/library/archaea
grep -o "ftp.*1.*$"  assembly_summary.txt |  awk '{print $1}' > archaea_ftps #lost a couple files but not going to worry about it
sed 's/ftp:/https:/g' archaea_ftps > archaea_https
cat archaea_https | sed 's:.*/::' > filenames
sed 's/$/_genomic/' filenames > filenames1
paste filenames1 archaea_https | column -s $'\t' -t > archaea_table

>grabs
while IFS=" " read -r a b
do
  echo "$b/$a.fna.gz" >> grabs
done < archaea_table

>get_genomes
for file in `cat grabs`; do echo "wget $file $NEW_DB/library/archaea/" >> get_genomes; done
ls5_launcher_creator.py -j get_genomes -n get_genomes -t 00:15:00 -N 3 -w 24 -a tagmap -e cbscott@utexas.edu -q normal

>unzipme
for file in GCF*.gz;
do echo "gunzip $file" >> unzipme; done
ls5_launcher_creator.py -j unzipme -n unzipme -t 00:15:00 -N 1 -w 24 -a tagmap -e cbscott@utexas.edu -q normal

#REMEMBER TO NAME JOBS WITH THE TAX NAME
#bacteria
cd $NEW_DB/library
mkdir bacteria
cd bacteria
#there's a ton of bacteria... want to give it more nodes and way more time.
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt $NEW_DB/library/bacteria
ls5_launcher_creator.py -j get_genomes -n bact_genomes -t 10:00:00 -N 15 -w 24 -a tagmap -e cbscott@utexas.edu -q normal

#invertebrates
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrates/assembly_summary.txt $DATABASE
ls5_launcher_creator.py -j get_genomes -n invert_genomes -t 02:30:00 -N 5 -w 24 -a tagmap -e cbscott@utexas.edu -q normal

ls5_launcher_creator.py -j get_genomes -n invert_genomes -t 14:00:00 -N 10 -w 20 -a tagmap -e cbscott@utexas.edu -q normal

#plants
ls5_launcher_creator.py -j get_genomes -n plant_genome -t 14:00:00 -N 8 -w 16 -a tagmap -e cbscott@utexas.edu -q normal

#protozoa
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/assembly_summary.txt $FILE_DATABASE
ls5_launcher_creator.py -j get_genomes -n get_genomes -t 00:30:00 -N 3 -w 24 -a tagmap -e cbscott@utexas.edu -q normal

>unzipme
for file in GCF*.gz;
do echo "gunzip $file" >> unzipme; done
ls5_launcher_creator.py -j unzipme -n unplant -t 00:10:00 -N 3 -w 24 -a tagmap -e cbscott@utexas.edu -q normal
ls5_launcher_creator.py -j unzipme -n uninvert -t 00:10:00 -N 1 -w 24 -a tagmap -e cbscott@utexas.edu -q development


ls5_launcher_creator.py -j totalref -n totalref -t 00:10:00 -N 1 -w 24 -a tagmap -e cbscott@utexas.edu -q development

#mammalian verts? Just humans, done
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz .


#Can we make this a shell script?
$ASSEMBLY_SUMMARY=$1
$DATABASE=$2
grep -o "ftp.*1.*$"  $ASSEMBLY_SUMMARY |  awk '{print $1}' > ftps #get ftps links
sed 's/ftp:/https:/g' ftps > https #make https links because ftp isn't working
cat https | sed 's:.*/::' > filenames #get just the filenames to finish the links
sed 's/$/_genomic/' filenames > filenames1 #add the necessary subfix to the links
paste filenames1 https | column -s $'\t' -t > table #create a master table
>grabs
while IFS=" " read -r a b #while the delimiter is space/tab, read in two columns
do
  echo "$b/$a.fna.gz" >> grabs #create the final https links
done < table
>get_genomes
for file in `cat grabs`; do echo "wget $file $DATABASE" >> get_genomes; done #wget these files

rm ftps
rm https
rm filenames1
rm filenames
rm table
rm grabs

#Keep wrangling with HOPS
#Hops jar file doesn't exist - get it from a different branch of rhuebler's git hub
wget https://github.com/rhuebler/HOPS/blob/external_update/AMPS/Resources/hops0.2.jar #replaced line in installation file with this
#Try an old HOPS from a past github tree
wget https://github.com/rhuebler/HOPS/blob/6734fe8776f690fb98030f9a72b44c2c0fe70e15/Resources/hops0.31.jar #doesn't work because need raw file
#download to computer and then scp
#Original: hops0.31.jar doesn't exist?
#https://raw.githubusercontent.com/rhuebler/HOPS/external/Resources/hops0.31.jar
#wget the tre and map Files
wget https://raw.githubusercontent.com/rhuebler/HOPS/external/Resources/ncbi.tre
wget https://raw.githubusercontent.com/rhuebler/HOPS/external/Resources/ncbi.map

ln -s $HOME/malt/malt-run $HOME/hops/malt-run #symbollic link to MALT Run #looks like there may be an issue
ln -s $SCRATCH/ancientDNA/malt-allmarine $HOME/hops/database #symbollic link to malt-built SILVA database

#let's test on SILVA db already constructed
echo "java -Xmx700G -jar $HOME/hops/hops.jar --input /scratch/06909/cbscott/ancientDNA/rawreads/S17466.Y1.E1.L1.fastq --output $SCRATCH/ancientDNA/HOPStest -m full" > hops_test
ls5_launcher_creator.py -j hops_test -n hops_test -t 00:20:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q development


#failed after MALT step.
#Try with config file?
echo "/home1/06909/cbscott/jdk-14.0.1/bin/java -Xmx700G -jar $HOME/hops/hops.jar --input /scratch/06909/cbscott/ancientDNA/rawreads/test.fastq --output $SCRATCH/ancientDNA/HOPStest -m full -c $SCRATCH/ancientDNA/configFile.txt" > hops_test_config
ls5_launcher_creator.py -j hops_test_config -n hops_test_config -t 01:00:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q normal #memory error here? this worked!!
#error

#Run MALT extract separately on existing rma file (memory error still occuring?)
java-new -jar  MaltExtract1.5.jar --help


#Missing list of interesting taxa?
#Attempt another MALT-build? This works; careful!!!!! acc2taxonomy has inconsistent notation in manual and help
#we need a master accession2taxid list because it seems like these are coming from more than one source
sed '1d' nucl_gb.accession2taxid > tmpfile; # POSIX
cat nucl_wgs.accession2taxid tmpfile > nucl_all.accession2taxid
#Also create a master genomic database (ACCIDENTALLY MOVED SOME FILES INSTEAD OF COPIED)

#decrease threads....
echo "malt-build --input $SCRATCH/ancientDNA/bact_refgenomes/plantless-malt/* -s DNA --index $SCRATCH/ancientDNA/malt-allmarine --acc2taxonomy $SCRATCH/ancientDNA/bact_refgenomes/NCBI_kraken/taxonomy/nucl_all.accession2taxid --threads 1 --verbose --step 5" > malt_noplant #still ran into a memory error
ls5_launcher_creator.py -j malt_noplant -n malt_noplant -t 01:00:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q development #memory error here? this worked!!

rm !(*.fa) #get rid of all not .fa
#Try to MALT build without plants... is it small enough?

#------------------Examine Damage patterns using PMDtools-----------------------
testbam=S17463.Y1.E1.L1.bam
echo "samtools view -h $testbam | pmdtools.0.60.py --threshold 1.5 --header | samtools view -Sb - > ${testbam/.bam/}.pmds1_5filter.bam" > pmds1_5filter
echo "samtools view -h $testbam | python pmdtools.0.60.py -p > pmdscores.txt" > pmdscores

ls5_launcher_creator.py -j pmds1_5filter -n pmds1_5filter -t 00:10:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q development

echo "samtools view $testbam | python pmdtools.0.60.py --platypus --requirebaseq 30 > PMD_temp.txt" > pmddist
ls5_launcher_creator.py -j pmddist -n pmddist -t 00:10:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q development


samtools view -h S17463.Y1.E1.L1.sam | pmdtools.0.60.py --deamination --header > PMD_temp.txt
samtools view -h S17463.Y1.E1.L1.pmds1_5filter.bam | pmdtools.0.60.py -p > postPMD_dist.txt #-h is to ensure the header is included in the output. pmdtools -p is for the pmdscore output

#copy files

for file in `cat MarRef`; do cp $file $SCRATCH/ancientDNA/bact_refgenomes/NCBI_kraken/library/MarRef/; done


>returnfiles
for line in filenames; do cut -c1-15 $line >> returnfiles; done

for file in `cat returnfiles`; do cp $SCRATCH/ancientDNA/bact_refgenomes/NCBI_kraken/library/total_reference/${file}* $SCRATCH/ancientDNA/bact_refgenomes/NCBI_kraken/library/protozoa/; done


echo "cp $SCRATCH/ancientDNA/bact_refgenomes/NCBI_kraken/library/MarRef/* $SCRATCH/ancientDNA/bact_refgenomes/plantless-malt" > copyfiles
ls5_launcher_creator.py -j copyfiles -n copyfiles -t 00:05:00 -N 1 -w 10 -a tagmap -e cbscott@utexas.edu -q development #memory error here?

echo "cp $SCRATCH/ancientDNA/bact_refgenomes/NCBI_kraken/library/protozoa/* $SCRATCH/ancientDNA/bact_refgenomes/plantless-malt" > copyfiles1
ls5_launcher_creator.py -j copyfiles1 -n copyfiles1 -t 00:10:00 -N 1 -w 10 -a tagmap -e cbscott@utexas.edu -q development #memory error here?


#Try KRAKEN2 classification?
echo "kraken2 --db $DBNAME /scratch/06909/cbscott/ancientDNA/rawreads/S17463.Y1.E1.L1.fastq" > kraken_classify
ls5_launcher_creator.py -j kraken_classify -n kraken_classify -t 00:10:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q development #memory error here?

#Try with a report flag?
echo "kraken2 --db $DBNAME /scratch/06909/cbscott/ancientDNA/rawreads/S17463.Y1.E1.L1.fastq --report" > kraken_classify_report
ls5_launcher_creator.py -j kraken_classify_report -n kraken_classify_report -t 00:10:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q development #memory error here?
#report didn't do anything?
ImportTaxomoy.pl -q 2 -t 3 YourKrakenOutputFile -o YourKronaReportFile #no matches :(

#try HOPS again with updated java?
#also, add SILVA to kraken DB and try again?
/home1/06909/cbscott/Krona/KronaTools/taxonomy/accession2taxid

#missing acesssion2taxids
echo "./updateAccessions.sh --only-build" > buildAcc
ls5_launcher_creator.py -j buildAcc -n buildAcc -t 01:00:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q normal

#need blast2rma help
echo "blast2rma -h" > helpme
ls5_launcher_creator.py -j helpme -n helpme -t 00:00:30 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q normal

#use MEGAN for viz!
#maybe need read file?
echo "blast2rma -i 464.blastn.gz -r $SCRATCH/ancientDNA/rawreads/464.fastq -f BlastText -o /scratch/06909/cbscott/ancientDNA/malt_tests -supp 0" > MEGANviz
ls5_launcher_creator.py -j MEGANviz -n MEGANviz -t 00:10:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q development


$WORK/Krona/KronaTools/scripts/ImportBLAST.pl -q 2 -t 3 test.blastn.tab -o testkrona.html -tax $WORK/Krona/KronaTools/taxonomy/

#Attempt to subdivide fastq file into 4 smaller FILES
split -l 50000000 S17464.Y1.E1.L1.fastq 464
ls 464a* > 464files

split -l 500000 S17464.Y1.E1.L1.fastq 464small
ls 464smalla* > 464smallfiles

#try just 464smallaa
echo "malt-run -i $SCRATCH/ancientDNA/rawreads/464smallaa -a $SCRATCH/ancientDNA/malt_tests --index $SCRATCH/ancientDNA/malt-allmarine/ --mode BlastN --format Text --alignmentType SemiGlobal" > small_file
ls5_launcher_creator.py -j small_file -n small_file -t 00:07:00 -a tagmap -e cbscott@utexas.edu -N 1 -w 1 -q development

#issues with blast2rma.... read file?
echo "blast2rma -i 464smallaa.blastn.gz -r $SCRATCH/ancientDNA/rawreads/464smallaa -f BlastText -o /scratch/06909/cbscott/ancientDNA/malt_tests -supp 0 -ms 50 -me 0.01 -mpi 95" > 464smallviz
ls5_launcher_creator.py -j 464smallviz -n 464smallviz -t 00:06:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q development


#try to malt extract 464smallaa ADJUST FILTERS?
echo "Citromicrobium
Erythrobacter flavus
Marinobacter hydrocarbonoclasticus" > bact_of_interest.txt
echo "$HOME/jdk-14.0.1/bin/java -jar  $HOME/hops/MaltExtract1.5.jar -i 464smallaa.rma6 -o $SCRATCH/ancientDNA/malt_tests/MALT_extract --resources /home1/06909/cbscott/hops --filter def_anc --taxa bact_of_interest.txt" > maltextract
ls5_launcher_creator.py -j maltextract -n maltextract -t 00:03:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q development


/home1/06909/cbscott/hops/postprocessing.AMPS.r -m def_anc -r $SCRATCH/ancientDNA/malt_tests/MALT_extract -t 1 -n $SCRATCH/ancientDNA/malt_tests/bact_of_interest.txt

#get algal genomes (all)
awk '{print $1}' assembly-summary > accessionNums
$HOME/datasets download assembly -i accessionNums -f algae_genomes.zip

#MUST BE IN PROPER DIRECTORY
find . -name '*.fna' -exec cp {} $SCRATCH/ancientDNA/bact_refgenomes/plantless-malt \;
#too many algae... try fewer. Delete all. Only take refseq and add in symbiodinium
find . -name '*.fna' > grabbed_files
sed 's|.*/||' grabbed_files > to_delete

for algae in `cat to_delete`; do rm /scratch/06909/cbscott/ancientDNA/bact_refgenomes/plantless-malt/${algae}; done

#subset algal genomes
awk '{print $1}' assembly_result-4.txt > accessionNums
$HOME/datasets download assembly -i accessionNums -f reduced_algae_genomes.zip
find . -name '*.fna' -exec cp {} $SCRATCH/ancientDNA/bact_refgenomes/plantless-malt \;

echo "malt-build --input $SCRATCH/ancientDNA/bact_refgenomes/plantless-malt/* -s DNA --index $SCRATCH/ancientDNA/malt-allmarine --acc2taxonomy $SCRATCH/ancientDNA/bact_refgenomes/NCBI_kraken/taxonomy/nucl_all.accession2taxid --threads 2 --step 15 --verbose" > malt_noplant
#No mem error after increasing step size... maybe 5 is too big. Try to decrease and refine in the future. AH! Forgot to add back *all* of the algae. If increasing step size seems to work, try adding back.
ls5_launcher_creator.py -j malt_noplant -n malt_noplant -t 00:30:00 -N 2 -w 1 -a tagmap -e cbscott@utexas.edu -q normal #memory error here? this worked!!

echo "malt-run -i $SCRATCH/ancientDNA/rawreads/subset/464smallaa -a $SCRATCH/ancientDNA/malt_tests --index $SCRATCH/ancientDNA/malt-allmarine/ --mode BlastN --format Text --alignmentType SemiGlobal" > algae_db_align
ls5_launcher_creator.py -j algae_db_align -n algae_db_align -t 00:30:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q development #memory error here? this worked!!

echo "blast2rma -i 464smallaa.blastn.gz -r $SCRATCH/ancientDNA/rawreads/subset/464smallaa -f BlastText -o /scratch/06909/cbscott/ancientDNA/malt_tests -supp 0" > MEGANviz
ls5_launcher_creator.py -j MEGANviz -n MEGANviz -t 00:10:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q development

echo "Symbiodinium
Symbiodinium kawagutii
Symbiodinium microadriaticum
Symbiodinium sp. clade A Y106
" > bact_of_interest.txt
#be sure to create a new file directory for the extraction- this will NOT overwrite old files.
echo "$HOME/jdk-14.0.1/bin/java -jar  $HOME/hops/MaltExtract1.5.jar -i 464smallaa.rma6 -o $SCRATCH/ancientDNA/malt_tests/extract_symbiodinium --resources /home1/06909/cbscott/hops --filter def_anc --taxa bact_of_interest.txt" > maltextract
ls5_launcher_creator.py -j maltextract -n maltextract -t 00:03:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q development

/home1/06909/cbscott/hops/postprocessing.AMPS.r -m def_anc -r $SCRATCH/ancientDNA/malt_tests/extract_symbiodinium -t 1 -n $SCRATCH/ancientDNA/malt_tests/bact_of_interest.txt
#getting weird errors with this for symbiodinium (too few ancient reads?)

#JUST ALGAE DATABASE? I.E., is there anything else but symbiodinium?
find . -name '*.fna' -exec cp {} $SCRATCH/ancientDNA/bact_refgenomes/algae_ref \; #takes 2-3 mins
echo "malt-build --input $SCRATCH/ancientDNA/bact_refgenomes/algal_ref -s DNA --index $SCRATCH/ancientDNA/malt-allalgae --acc2taxonomy $SCRATCH/ancientDNA/bact_refgenomes/NCBI_kraken/taxonomy/nucl_all.accession2taxid --threads 1 --step 1 --verbose" > malt_allalgae
ls5_launcher_creator.py -j malt_allalgae -n malt_allalgae -t 00:25:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q normal

echo "malt-run -i $SCRATCH/ancientDNA/rawreads/subset/464smallaa -a $SCRATCH/ancientDNA/malt_tests --index $SCRATCH/ancientDNA/malt-allalgae/ --mode BlastN --format Text --alignmentType SemiGlobal" > small_file
ls5_launcher_creator.py -j small_file -n small_file -t 00:07:00 -a tagmap -e cbscott@utexas.edu -N 1 -w 1 -q development

echo "blast2rma -i 464smallaa.blastn.gz -r $SCRATCH/ancientDNA/rawreads/subset/464smallaa -f BlastText -o /scratch/06909/cbscott/ancientDNA/malt_tests -supp 0" > MEGANviz
ls5_launcher_creator.py -j MEGANviz -n MEGANviz -t 00:05:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q development
#just algae does a bad job! does not map to something that makes sense.

#Can we use SILVA? I feel like I'm going full circle, but worth a try.
echo "malt-build --input /scratch/06909/cbscott/ancientDNA/bact_refgenomes/silva_dna/* -s DNA --index $SCRATCH/ancientDNA/malt-SILVA --acc2taxonomy $SCRATCH/ancientDNA/bact_refgenomes/NCBI_kraken/taxonomy/nucl_all.accession2taxid --step 1 --verbose" > malt_SILVA
ls5_launcher_creator.py -j malt_SILVA -n malt_SILVA -t 01:00:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q normal

echo "malt-run -i $SCRATCH/ancientDNA/rawreads/subset/464smallaa -a $SCRATCH/ancientDNA/malt_tests --index $SCRATCH/ancientDNA/malt-SILVA/ --mode BlastN --format Text --alignmentType SemiGlobal" > small_file
ls5_launcher_creator.py -j small_file -n small_file -t 00:07:00 -a tagmap -e cbscott@utexas.edu -N 1 -w 1 -q development

echo "blast2rma -i 464smallaa.blastn.gz -r $SCRATCH/ancientDNA/rawreads/subset/464smallaa -f BlastText -o /scratch/06909/cbscott/ancientDNA/malt_tests -supp 0" > MEGANviz
ls5_launcher_creator.py -j MEGANviz -n MEGANviz -t 00:05:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q development

#How well does mapping do just to MarRef? JUST A MARREF DATABASE. Can we get away with merging algae and marref? This is very large...
echo "malt-build --input /scratch/06909/cbscott/ancientDNA/bact_refgenomes/NCBI_kraken/library/MarRef/* -s DNA --index $SCRATCH/ancientDNA/malt-MarRef --acc2taxonomy $SCRATCH/ancientDNA/bact_refgenomes/NCBI_kraken/taxonomy/nucl_all.accession2taxid --threads 1 --step 1 --verbose" > malt_MarRef
ls5_launcher_creator.py -j malt_MarRef -n malt_MarRef -t 00:30:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q normal

#What about SILVA?

#deal with database issues: add algae, remove bacteria that nothing mapped to.... try CENTRIFUGE? Uses the entire ncbi nt database. Then pipe into malt_extract
#end of day 06/23:
#Re-map files
#Grab nt database, and build?
wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz
echo "gunzip nt.gz && mv -v nt nt.fa" > unzip_ntdb #too big, needs to be a job :( 15 minutes not enough! ah!
ls5_launcher_creator.py -j unzip_ntdb -n unzip_ntdb -t 02:00:00 -N 1 -w 1 -a tagmap -e cbscott@utexas.edu -q normal

# Get mapping file
#wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz original from manual
wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz #try this, use nucleotide acessions2taxid instead of gi_taxid_nucl map....
gunzip -c taxdump.tar.gz | sed 's/^/gi|/' > gi_taxid_nucl.map

# build index using 16 cores and a small bucket size, which will require less memory (swap out mapping file)
echo "centrifuge-build -p 16 --bmax 1342177280 --conversion-table all_nucl.accession2taxid \
                 --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
                 nt.fa nt" > test_centri_build
ls5_launcher_creator.py -j test_centri_build -n test_centri_build -t 01:00:00 -N 2 -w 8 -a tagmap -e cbscott@utexas.edu -q normal


#I think centrifuge will output a .sam file, that can then be converted to an rma6 file and converted for use in MEGAN and HOPS
>maps_czoox
for file in *.fastq; do
echo "bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x $GENOME_FASTA -U $file -S ${file/.fastq/}.sam && \
samtools sort -O bam -o ${file/.fastq/}.coralzoox.sorted.bam ${file/.fastq/}.sam && samtools index -c ${file/.fastq/}.coralzoox.sorted.bam " >> maps_czoox;
done
ls5_launcher_creator.py -j maps_czoox -n maps_czoox -t 2:00:00 -N 2 -w 24 -a tagmap -e cbscott@utexas.edu -q normal

#break it into smaller jobs. That way they sit in queue for less time?
export GENOME_FASTA=$WORK/ancientDNAdb/Apalm_sym_euk.fasta
>maps_czoox_1
for file in *.fastq; do
echo "bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x $GENOME_FASTA -U $file -S ${file/.fastq/}.sam">>maps_czoox_1; done
ls5_launcher_creator.py -j maps_czoox_1 -n maps_czoox_1 -t 02:00:00 -N 2 -w 24 -a tagmap -e cbscott@utexas.edu -q normal

>maps_czoox_2
for file in *.fastq; do
echo "samtools sort -O bam -o ${file/.fastq/}.coralzoox.sorted.bam ${file/.fastq/}.sam" >> maps_czoox_2; done
ls5_launcher_creator.py -j maps_czoox_2 -n maps_czoox_2 -t 00:30:00 -N 1 -w 24 -a tagmap -e cbscott@utexas.edu -q development
