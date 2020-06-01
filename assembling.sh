#This script records bash commands used for fungal genome assembly

##acquire information about all subreads and filtered subreads in the dataset
ls -lah ../../shared_bioinformatics_master_projects/agaricalesGenomes/b2016040/INBOX/*/*subreads/;
##create corresponding repositories in the fastqcReports directory
ls ../../shared_bioinformatics_master_projects/agaricalesGenomes/b2016040/INBOX/*/*subreads/*|cut -d / -f 7|while read fungi_name;do mkdir ./fastqcReports/$fungi_name;done;
##Do fastqc analysis for the 4 given fastq.gz files
ls fastqcReports/*|cut -d / -f 2|sed s/://g|while read name;do nohup fastqc --threads 20 -o ./fastqcReports/$name/ --extract -f fastq -c ../../shared_bioinformatics_master_projects/agaricalesGenomes/b2016040/INBOX/$name/*reads/*.gz;done;
#Do multiqc analysis from results of fastqcReports
multiqc ./fastqcReports/ -o multiQCreport/;
##From the supplementary material it can be learned that the genome sequencing and RNA sequencing were done by PacBio, so '--pacbio' option should be used
##create corresponding repositories in the spadesAssembly directory
ls ../../shared_bioinformatics_master_projects/agaricalesGenomes/b2016040/INBOX/*/*subreads/*|cut -d / -f 7|while read fungi_name;do mkdir ./spadesAssembly/$fungi_name;done;

##Create .gitignore
ls -d */|sed "s/\///g"|while read name;do echo $name >> .gitignore;done;

##Long PacBio RSII reads assembly by wtdbg2 to create concensus genomes for further genome polishing by arrow 
nohup ./wtdbg2 -x rs -g ? -i ../pb_279_filtered_subreads.fastq.gz -t 20 -fo ../wtdbg2Assembly/pb_279_default/pb_279_default;
nohup ./wtdbg2 -x rs -g ? -i ../pb_320-2_filtered_subreads.fastq.gz -t 20 -fo ../wtdbg2Assembly/pb_320-2_default/pb_320-2_default;
# derive raw consensus
ls ../wtdbg2Assembly/*/*.ctg.lay.gz|while read path;do dir=$(echo $path|cut -d / -f 1-3);name=$(echo $path|cut -d / -f 4|cut -d . -f 1);nohup ./wtpoa-cns -t 20 -i $path -fo $dir/$name.raw.fa;done;

##In order to align original .fastq.gz files of long reads into .BAM files for genome polishing with arrow, it's necessary to use blasr(bowtie2 is used for short reads, inappropriate here)
ls ../unzipped/*.fastq|while read subreads;do name=$(echo $subreads|cut -d / -f 3|cut -d _ -f 1,2);ref=$(echo ./$name*/$name*.raw.fa);dir=$(echo $ref|cut -d / -f 1,2);nohup blasr $subreads $ref --bam --out $dir/blasrBAM/$name.bam;done  #however, blasr doesn't work here, 'your subreads aren't PacBio reads'
# As an attempt to solve the problem above, convert .fastq.gz files to .fasta subread files
nohup ls *_filtered_subreads.fastq.gz|while read gz;do name=$(echo $gz|cut -d . -f 1);less $gz|grep -A 1 ^"@"|sed -r 's/^@/>/'> ${name}.fasta1;echo $name;done
ls *.fasta1|while read fasta1;do name=$(echo $fasta1|cut -d . -f 1);grep -v "--" < $fasta1 > ${name}.fasta;done

##Long read assembly by Miniasm
# Find overlaps by all-vs-all self-mappings via minimap2
nohup minimap2 -x ava-pb -t 15 pb_279_filtered_subreads.fastq.gz pb_279_filtered_subreads.fastq.gz|gzip -1 >miniasmAssembly/pb_279_default.paf.gz
nohup minimap2 -x ava-pb -t 15 pb_320-2_filtered_subreads.fastq.gz pb_320-2_filtered_subreads.fastq.gz|gzip -1 >miniasmAssembly/pb_320-2_default.paf.gz
# Layout using Miniasm
ls miniasmAssembly/*|while read paf;do name=$(echo $paf|cut -d / -f 2|cut -d . -f 1|cut -d _ -f 1,2);nohup miniasm -f $name*.fastq.gz $paf > ./miniasmAssembly/$name.gfa;done 
# Convert .gfa to .fasta
awk '/^S/{print ">"$2"\n"$3}' pb_279.gfa > pb_279-unpolished.fasta
awk '/^S/{print ">"$2"\n"$3}' pb_320-2.gfa > pb_320-2_unpolished.fasta

# Assembly polishing by minipolish(the output files aren't actually in .gfa format, eerie here)
nohup minipolish -t 20 --pacbio pb_279_filtered_subreads.fastq.gz miniasmAssembly/pb_279.gfa > miniasmAssembly/pb_279_polished.gfa
nohup minipolish -t 20 --pacbio pb_320-2_filtered_subreads.fastq.gz miniasmAssembly/pb_320-2.gfa > miniasmAssembly/pb_320-2_polished.gfa
# Convert .gfa to .fasta
awk '/^S/{print ">"$2"\n"$3}' pb_279_polished.gfa > pb_279-polished.fasta
awk '/^S/{print ">"$2"\n"$3}' pb_320-2_polished.gfa > pb_320-2_polished.fasta

## Genome assembling and polishing by raven with default setting
ls *.fastq.gz|while read fastq;do name=$(echo $fastq|cut -d . -f 1);nohup raven --graphical-fragment-assembly ravenAssembly/default/$name.gfa  -t 20  $fastq;echo $name;done 
# change the alignment parameters
ls *.fastq.gz|while read fastq;do name=$(echo $fastq|cut -d . -f 1);nohup raven -m 2 -n -6 -g -4 --graphical-fragment-assembly ravenAssembly/par1/${name}p2m2n-6g-4.gfa -t 20 $fastq;echo $fastq;done
ls *.fastq.gz|while read fastq;do name=$(echo $fastq|cut -d . -f 1);nohup raven -p 3 -m 3 -n -5 -g -4 --graphical-fragment-assembly ravenAssembly/par2/${name}p3m3n-5g-4.gfa -t 20 $fastq;echo $fastq;done;
ls *.fastq.gz|while read fastq;do name=$(echo $fastq|cut -d . -f 1);nohup raven -p 3 -m 2 -n -6 -g -4 --graphical-fragment-assembly ravenAssembly/par3/${name}p3m2n-6g-4.gfa -t 20 $fastq;echo $fastq;done;
ls *.fastq.gz|while read fastq;do name=$(echo $fastq|cut -d . -f 1);nohup raven -p 1 -m 2 -n -6 -g -4 --graphical-fragment-assembly ravenAssembly/par4/${name}p1m2n-6g-4.gfa -t 20 $fastq;echo $fastq;done;
ls *.fastq.gz|while read fastq;do name=$(echo $fastq|cut -d . -f 1);nohup raven -p 1 -m 3 -n -5 -g -4 --graphical-fragment-assembly ravenAssembly/par5/${name}p1m3n-5g-4.gfa -t 20 $fastq;echo $fastq;done;
# Convert .gfa files to .fasta files
awk '/^S/{print ">"$2"\n"$3}' pb_279_filtered_subreads.gfa >pb_279.fasta;
awk '/^S/{print ">"$2"\n"$3}' pb_320-2_filtered_subreads.gfa >pb_320-2.fasta; 

# Generate Bandage plots for .gfa assemblies
ls *.gfa|while read gfa;do Bandage image $gfa $gfa.jpg;done

## canu(subreads correction & trimming & genome assembly)
nohup canu -p pb_279 -d canuAssembly/default/pb_279/ genomeSize=74m errorRate=0.3 gnuplotTested=true -pacbio-raw pb_279_filtered_subreads.fastq.gz 
nohup canu -p pb_320-2 -d canuAssembly/default/pb_320-2/ genomeSize=137m errorRate=0.3 gnuplotTested=true -pacbio-raw pb_320-2_filtered_subreads.fastq.gz 

## repeat graph genome assembly & polishing by flye
nohup flye --pacbio-raw pb_279_filtered_subreads.fastq.gz --genome-size 74m --out-dir flyeAssembly/default/pb_279/ --threads 20
nohup flye --pacbio-raw pb_320-2_filtered_subreads.fastq.gz --genome-size 137m --out-dir flyeAssembly/default/pb_320-2/ --threads 20

## Quast analysis for genome assemblies
## As we know, the minimum threshold of a contig length for analysis in Quast's default setting is 500. However, this might not be the case here. Distribution of contig lengths of the previously generated genome assemblies has to be analyzed here using awk.
# Check raven & miniasm assemblies ./default
ls *.fasta|while read file;do for thres in {100..500..100};do awk -v file="$file" -v i="$thres" '/^[^>]/ {sum++;if(length($0)>=i){qualified++}} END{print qualified/sum}' $file;done;done #All output values are 1, 500 can be unchanged in QUAST analysis
# concatenate original genome assemblies to create reference genomes 
cat *.fa > ../../../../../jiawei_zhao/LUfungiProject/OriginalAssemblies/pb_320-2_Mysco.fasta
cat *.fa > ../../../../../jiawei_zhao/LUfungiProject/OriginalAssemblies/pb_279_Leuge.fasta
# Quast analysis
nohup quast.py -o pb_320-2_quast  pb_320-2_polished.fasta -r ../../OriginalAssemblies/pb_320-2_Mysco.fasta -t 20
nohup quast.py -o pb_279_quast/ pb_279-polished.fasta -r ../../OriginalAssemblies/pb_279_Leuge.fasta -t 20
# find out all directories which store the results of Quast analysis
find -type d -name "*quast*"

## Kmerfreq analysis for files of subreads
nohup ./kmerfreq -f 1 -p pb_279 -t 30 subreads1.lib;
nohup ./kmerfreq -f 1 -t 30 -p pb_320-2 subreads2.lib; 

##Assembly by pb-assembly
#Convert .fastq.gz files to .fasta files as required
nohup less pb_279_filtered_subreads.fastq.gz|paste - - - -|sed 's/^@/>/'|awk '{print $1"\n"$2}' > pb_279_filtered_subreads.fasta;
nohup less pb_320-2_filtered_subreads.fastq.gz|paste - - - -|sed 's/^@/>/'|awk '{print $1"\n"$2}' > pb_320-2_filtered_subreads.fasta;
# Create symbolic links for the directories which contain .bax.h5 files
ln -s ../../shared_bioinformatics_master_projects/agaricalesGenomes/b2016040/INBOX/pb_279/rawdata/run1 pb_279_raw
ln -s ../../shared_bioinformatics_master_projects/agaricalesGenomes/b2016040/INBOX/pb_320-2/raw/run1 pb_320-2_raw
#Create an input list for .bax.h5 files
ls pb_279_raw/*/Analysis_Results/*.bax.h5 >> pb_279_baxh5.txt
ls pb_320-2_raw/*/Analysis_Results/*.bax.h5 >> pb_320-2_baxh5.txt
for i in {3..24..3};do sed -n "$(echo "$i-2"|bc),${i}p" < pb_279_baxh5.txt >> bamfiles/pb_279_list/pb_279_bax_list$(echo "$i/3"|bc).txt;done
for i in {3..24..3};do sed -n "$(echo "$i-2"|bc),${i}p" < pb_320-2_baxh5.txt >> bamfiles/pb_320-2_list/pb_320-2_bax_list$(echo "$i/3"|bc).txt;done
#Use bax2bam to create .bam subread files which are necessary for input of pb-assembly
for i in {1..8..1};do nohup bax2bam -f bamfiles/pb_279_list/pb_279_bax_list${i}.txt -o bamfiles/pb_279/${i}subreads --subread;done
for i in {1..8..1};do nohup bax2bam -f bamfiles/pb_320-2_list/pb_320-2_bax_list${i}.txt -o bamfiles/pb_320-2/${i}subreads --subread;done 
#Create .fofn files which contain paths of .fasta and .bam files for falcon & falcon-unzip
cd pb-assembly
ls ../bamfiles/pb_279/*.subreads.bam|while read bam;do echo $bam >> pb_279_bam.fofn;done
ls ../bamfiles/pb_320-2/*.subreads.bam|while read bam;do echo $bam >> pb_320-2_bam.fofn;done
echo "../pb_279_filtered_subreads.fasta" >> pb_279_fa.fofn 
echo "../pb_320-2_filtered_subreads.fasta" >> pb_320-2_fa.fofn 
#Download reference .cfg files and then modify it to fit specific requirements & download evaluation scripts
wget https://pb-falcon.readthedocs.io/en/latest/_downloads/fc_run_fungal.cfg
git init 
git remote add cfg https://github.com/PacificBiosciences/pb-assembly/
git config core.sparsecheckout true
echo "cfgs" >> .git/info/sparse-checkout
echo "scripts" >> .git/info/sparse-checkout #only these two repositories ("cfgs" and "scripts") would be added into the checkout list of pulling
git pull cfg master
cp cfgs/fc_unzip.cfg mycfgs/
cp cfgs/fc_phase.cfg mycfgs/
#Start fc_run to do 0.pre-assembly 1.pread overlapping 2.contig assembly 
nohup fc_run ../mycfgs/fc_pb_279_v1.cfg &> run0.log &  #the first &> means to output both stderr and stdout to run0.log, & at the end of this command means running the command on background
nohup fc_run ../mycfgs/fc_pb_320-2_v1.cfg &> run0.log & 
#Evaluate Assembly Performance
python ../scripts/get_asm_stats.py 2-asm-falcon/p_ctg.fa
#Run FALCON-unzip to do 3.unzip 4.polish
mv "all.log" "all0.log"
nohup fc_unzip.py ../mycfgs/fc_unzip_pb_279.cfg &> run1.std &
nohup fc_unzip.py ../mycfgs/fc_unzip_pb_320-2.cfg &> run1.std &

