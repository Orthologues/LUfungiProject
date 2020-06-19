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
#Create its corresponding env
conda create -n denovo
conda activate denovo
conda install -c bioconda pb-assembly
conda install -c bioconda samtools=1.9 --force-reinstall #make it compatible with openssl 1.1.1
#Convert .fastq.gz files to .fasta files as required
nohup less pb_279_filtered_subreads.fastq.gz|paste - - - -|sed 's/^@/>/'|awk '{print $1"\n"$2}' > pb_279_filtered_subreads.fasta;
nohup less pb_320-2_filtered_subreads.fastq.gz|paste - - - -|sed 's/^@/>/'|awk '{print $1"\n"$2}' > pb_320-2_filtered_subreads.fasta;
#Calculate the sequence length distribution of .fasta files
nohup ls *.fasta|while read file;do nohup awk -v file="$file"  '/^[^>]/ {sum++;for(i=2000;i<=20000;i+=1000){if(length($0)>=i){qualified[(i-2000)/1000]++}}} END{print sum >> "lenstat.txt";for(i=0;i<=18;i++){printf("%d - %.2f\n",i*1000+2000,qualified[i]/sum*100) >> "lenstat.txt"}}' $file;done &
# Create symbolic links for the directories which contain .bax.h5 files
ln -s ../../shared_bioinformatics_master_projects/agaricalesGenomes/b2016040/INBOX/pb_279/rawdata/run1 pb_279_raw
ln -s ../../shared_bioinformatics_master_projects/agaricalesGenomes/b2016040/INBOX/pb_320-2/raw/run1 pb_320-2_raw
#Create an input list for .bax.h5 files
ls pb_279_raw/*/Analysis_Results/*.bax.h5 >> pb_279_baxh5.txt
ls pb_320-2_raw/*/Analysis_Results/*.bax.h5 >> pb_320-2_baxh5.txt
for i in {3..24..3};do sed -n "$(echo "$i-2"|bc),${i}p" < pb_279_baxh5.txt >> bamfiles/pb_279_list/pb_279_bax_list$(echo "$i/3"|bc).txt;done
for i in {3..24..3};do sed -n "$(echo "$i-2"|bc),${i}p" < pb_320-2_baxh5.txt >> bamfiles/pb_320-2_list/pb_320-2_bax_list$(echo "$i/3"|bc).txt;done
#Use bax2bam to create .bam subread files which are necessary for input of pb-assembly
for i in {1..8..1};do nohup bax2bam -f bamfiles/pb_279_list/pb_279_bax_list${i}.txt -o bamfiles/pb_279/${i} --subread;done
for i in {1..8..1};do nohup bax2bam -f bamfiles/pb_320-2_list/pb_320-2_bax_list${i}.txt -o bamfiles/pb_320-2/${i} --subread;done 
#Create .fofn files which contain paths of .fasta and .bam files for falcon & falcon-unzip
cd pb-assembly
ls ../bamfiles/pb_279/*subreads.bam|while read bam;do echo $bam >> pb_279_bam.fofn;done
ls ../bamfiles/pb_320-2/*subreads.bam|while read bam;do echo $bam >> pb_320-2_bam.fofn;done
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
nohup python ../scripts/get_asm_stats.py 2-asm-falcon/p_ctg.fa
#Create symbolic links of original unpolished contig assemblies in my repo
ln -s ../../../../shared_bioinformatics_master_projects/agaricalesGenomes/genome_assemblies/pb_320-2_Mysco/ ../../OriginalAssemblies/
ln -s ../../../../shared_bioinformatics_master_projects/agaricalesGenomes/genome_assemblies/pb_279_Leuge/ ../../OriginalAssemblies/
#Run FALCON-unzip to do 3.unzip 4.polish
mv "all.log" "all0.log"
nohup fc_unzip.py ../mycfgs/fc_unzip_pb_279.cfg &> run1.std &
nohup fc_unzip.py ../mycfgs/fc_unzip_pb_320-2.cfg &> run1.std & 
wait
#However, though '3-unzip' step was successful, '4-polish' failed because of the same unsolved issue as what the link https://github.com/PacificBiosciences/FALCON_unzip/issues/159 addresses. I decided to concatenate the output primary-contig & the associate-contig output file of 2-asm and the primary-contig & the haplotype-contig output file of '3-unzip' step in order to compare them to the original assembly each. Thus, blasr and arrow will be used later to polish the concatenated output file of '3-unzip' step.
cat 2-asm-falcon/p_ctg.fasta 2-asm-falcon/a_ctg.fasta > pb_320-2_falcon_step2_v1.fasta
cat 3-unzip/all_p_ctg.fasta 3-unzip/all_h_ctg.fasta >pb_320-2_falcon_step3_v1.fasta
cat 2-asm-falcon/p_ctg.fasta 2-asm-falcon/a_ctg.fasta > pb_279_falcon_step2_v1.fasta
cat 3-unzip/all_p_ctg.fasta 3-unzip/all_h_ctg.fasta >pb_279_falcon_step3_v1.fasta
cp *.fasta ../../../../shared_bioinformatics_master_projects/agaricalesGenomes/jiawei_zhao_assemblies/pb-assembly/
cp *.fasta ../../../../shared_bioinformatics_master_projects/agaricalesGenomes/jiawei_zhao_assemblies/pb-assembly/
#do QUAST analysis for comparison 
nohup quast.py -o step2_v1_quast/  pb_279_falcon_step2_v1.fasta -r ../../OriginalAssemblies/pb_279_Leuge.fasta -t 20 &
nohup quast.py -o step3_v1_quast/  pb_279_falcon_step3_v1.fasta -r ../../OriginalAssemblies/pb_279_Leuge.fasta -t 20 &
nohup quast.py -o step2_v1_quast/ pb_320-2_falcon_step2_v1.fasta -r ../../OriginalAssemblies/pb_320-2_Mysco.fasta -t 20 &
nohup quast.py -o step3_v1_quast/ pb_320-2_falcon_step3_v1.fasta -r ../../OriginalAssemblies/pb_320-2_Mysco.fasta -t 20 &
wait
find -maxdepth 3 -name "*.pdf"|while read pdf;do cp $pdf ../../../shared_bioinformatics_master_projects/agaricalesGenomes/jiawei_zhao_assemblies/pb-assembly/;done
#run blasr to generate aligned .bam files
cd ~/LUfungiProject/pb-assembly/pb_279_v1
ls ../../bamfiles/pb_279/*subreads.bam|while read bam;do echo $bam >> pb_279_bam.fofn;done
nohup blasr pb_279_bam.fofn pb_279_falcon_step3_v1.fasta --bam --out pb_279_v1_aligned.bam --nproc 20 &
cd ~/LUfungiProject/pb-assembly/pb_320-2_v1
ls ../../bamfiles/pb_320-2/*subreads.bam|while read bam;do echo $bam >> pb_320-2_bam.fofn;done
nohup blasr pb_320-2_bam.fofn pb_320-2_falcon_step3_v1.fasta --bam --out pb_320-2_v1_aligned.bam --nproc 20 & 
#blasr is outdated and doesn't work here try pbmm2 instead
#In order to make the output files compatible with genomicconsensus, concatenate separate .bam files to one .bam file for each fungus
cd bamfiles
nohup cat pb_279/*subreads.bam > pb_279_subreads.bam &
nohup cat pb_320-2/*subreads.bam > pb_320-2_subreads.bam &
#run pbmm2 to generate aligned .bam files
nohup pbmm2 align --preset SUBREAD pb_279_falcon_step3_v1.fasta ../../bamfiles/pb_279_subreads.bam pb_279_v1_step3_aligned.bam &
nohup pbmm2 align --preset SUBREAD pb_320-2_falcon_step3_v1.fasta ../../bamfiles/pb_320-2_subreads.bam pb_320-2_v1_step3_aligned.bam &
#pbmm2 doesn't work as well. Error message is 'terminate called without an active exception'. 

#Try to recreate 'denovo' env with py2.7
conda env remove -n denovo
conda create -n denovo python=2.7
conda activate denovo
conda install -c bioconda/label/cf201901 nim-falcon #https://github.com/PacificBiosciences/pbbioconda/issues/265
conda install -c bioconda pb-assembly
conda install -c bioconda samtools=1.9 --force-reinstall #https://github.com/bioconda/bioconda-recipes/issues/13958
#Change parameters to v2 and run falcon-assembler again
nohup fc_run ../mycfgs/fc_pb_279_v2.cfg &> run0.log &
nohup fc_run ../mycfgs/fc_pb_320-2_v2.cfg &> run0.log &
wait
nohup fc_unzip.py ../mycfgs/fc_unzip_pb_279.cfg &> run1.std &
nohup fc_unzip.py ../mycfgs/fc_unzip_pb_320-2.cfg &> run1.std & 
wait
#However, though '3-unzip' step was successful, '4-polish' failed because of the same unsolved issue as what the link https://github.com/PacificBiosciences/FALCON_unzip/issues/159 addresses. I decided to concatenate the output primary-contig & the associate-contig output file of 2-asm and the primary-contig & the haplotype-contig output file of '3-unzip' step in order to compare them to the original assembly each. Thus, blasr and arrow will be used later to polish the concatenated output file of '3-unzip' step.
# In order to save disk space, extract only .fa files out and delete the intermediate folders
mv 2-asm-falcon/*.fa .
mv 3-unzip/*.fa .
find -maxdepth 1|grep \.\/[0-9].|while read dir;do rm -rf $dir;done
#Generate .bam alignment files
nohup blasr ../../bamfiles/pb_279_subreads.bam pb_279_falcon_step3_v2.fasta --bam --out pb_279_step3_v2_aligned.bam --nproc 20 &
nohup blasr ../../bamfiles/pb_320-2_subreads.bam pb_320-2_falcon_step3_v2.fasta --bam --out pb_320-2_step3_v2_aligned.bam --nproc 20 & 
# all these blasr commands terminated with only output files in .bam.tmp format. From msg in nohup.out, it seems that the unaligned .bam files which I generated above have problems. Nonetheless, I changed their suffix from .bam.tmp to .bam and went forward as an attempt.
nohup arrow pb_279_step3_v2_aligned.bam -r pb_279_falcon_step3_v2.fasta -o pb_279_step3_v2_polished.fasta -o pb_279_step3_v2_polished.fastq -j 20 -pdb --diploid &
nohup arrow pb_320-2_step3_v2_aligned.bam -r pb_320-2_falcon_step3_v2.fasta -o pb_320-2_step3_v2_polished.fasta -o pb_320-2_step3_v2_polished.fastq -j 20 -pdb --diploid &
# However, 'no BGZF EOF marker; file may be truncated' is given here from genomicconsensus. In order to solve this issue, I used the script from https://github.com/peterjc/picobio/blob/master/sambam/bgzf_add_eof.py to process .bam files
./bgzf_add_eof.py pb_279_step3_v2_aligned.bam
./bgzf_add_eof.py pb_320-2_step3_v2_aligned.bam 
# Try genomicconsensus again with the revised .bam files
nohup arrow pb_279_step3_v2_aligned.bam -r pb_279_falcon_step3_v2.fasta -o pb_279_step3_v2_polished.fasta -o pb_279_step3_v2_polished.fastq -j 20 -pdb --diploid &
nohup arrow pb_320-2_step3_v2_aligned.bam -r pb_320-2_falcon_step3_v2.fasta -o pb_320-2_step3_v2_polished.fasta -o pb_320-2_step3_v2_polished.fastq -j 20 -pdb --diploid &
# An eerie error occured though: [ERROR] Genomic Consensus only works with cmp.h5 files and BAM files with accompanying .pbi files

# Try to generate .bam files from .bax.h5 files again
cd ~/LUfungiProject/
nohup sh -c 'for i in {1..8..1};do nohup bax2bam -f bamfiles/pb_279_list/pb_279_bax_list${i}.txt -o bamfiles/pb_279/${i} --subread --allowUnrecognizedChemistryTriple;done' &
nohup sh -c 'for i in {1..8..1};do nohup bax2bam -f bamfiles/pb_320-2_list/pb_320-2_bax_list${i}.txt -o bamfiles/pb_320-2/${i} --subread --allowUnrecognizedChemistryTriple;done' &
wait
cd ~/LUfungiProject/pb-assembly/
rm pb_279_bam.fofn
rm pb_320-2_bam.fofn
# Try assembly polishing, major steps are shown in https://www.biostars.org/p/273447/
cd ~/LUfungiProject/pb-assembly/pb_279_v2
ls ../../bamfiles/pb_279/*subreads.bam|while read bam;do echo $bam >> pb_279_bam.fofn;done
nohup pbmm2 align pb_279_falcon_step3_v2.fasta pb_279_bam.fofn pb_279_v2_aligned.bam --sort -j 8 -J 8 -m 32G --preset SUBREAD & #pbmm2 succeeds this time
nohup samtools faidx pb_279_falcon_step3_v2.fasta -o pb_279_falcon_step3_v2.fasta.fai &
nohup pbindex pb_279_v2_aligned.bam &
nohup arrow pb_279_v2_aligned.bam -r pb_279_falcon_step3_v2.fasta -o pb_279_step3_v2_polished.fastq -j 20 --diploid & wait
nohup cat pb_279_step3_v2_polished.fastq|paste - - - -|sed 's/^@/>/'|awk '{print $1"\n"$2}' > pb_279_step3_v2_polished.fasta & wait
cd ~/LUfungiProject/pb-assembly/pb_320-2_v2
ls ../../bamfiles/pb_320-2/*subreads.bam|while read bam;do echo $bam >> pb_320-2_bam.fofn;done
nohup pbmm2 align pb_320-2_falcon_step3_v2.fasta pb_320-2_bam.fofn pb_320-2_v2_aligned.bam --sort -j 8 -J 8 -m 32G --preset SUBREAD &
nohup samtools faidx pb_320-2_falcon_step3_v2.fasta -o pb_320-2_falcon_step3_v2.fasta.fai &
nohup pbindex pb_320-2_v2_aligned.bam &
nohup arrow pb_320-2_v2_aligned.bam -r pb_320-2_falcon_step3_v2.fasta -o pb_320-2_step3_v2_polished.fastq -j 20 --diploid &
nohup cat pb_320-2_step3_v2_polished.fastq|paste - - - -|sed 's/^@/>/'|awk '{print $1"\n"$2}' > pb_320-2_step3_v2_polished.fasta &

#create a more automated analysis workflow for assembly polishing
touch ~/LUfungiProject/countDone.txt
cd ~/LUfungiProject/pb-assembly/
versions1=$(find -maxdepth 1 -name "pb_279_v*"|wc -l);
versions2=$(find -maxdepth 1 -name "pb_320-2_v*"|wc -l);
for ((i=1;i<=$versions1;i++))
do
( cd ~/LUfungiProject/pb-assembly/pb_279_v${i}
  cat 2-asm-falcon/p_ctg.fa 2-asm-falcon/a_ctg.fa > pb_279_falcon_step2_v${i}.fasta
  cat 3-unzip/all_p_ctg.fa 3-unzip/all_h_ctg.fa >pb_279_falcon_step3_v${i}.fasta
  mv 2-asm-falcon/*.fa .
  mv 3-unzip/*.fa .
  find -maxdepth 1|grep \.\/[0-9].|while read dir;do rm -rf $dir;done
  ls ../../bamfiles/pb_279/*subreads.bam|while read bam;do echo $bam >> pb_279_bam.fofn;done
  nohup pbmm2 align pb_279_falcon_step3_v${i}.fasta pb_279_bam.fofn pb_279_v${i}_aligned.bam --sort -j 8 -J 8 -m 32G --preset SUBREAD & 
  nohup samtools faidx pb_279_falcon_step3_v${i}.fasta -o pb_279_falcon_step3_v${i}.fasta.fai &
  wait
  nohup pbindex pb_279_v${i}_aligned.bam &
  wait
  nohup arrow pb_279_v${i}_aligned.bam -r pb_279_falcon_step3_v${i}.fasta -o pb_279_step3_v${i}_polished.fastq -j 20 --diploid &
  wait
  nohup cat pb_279_step3_v${i}_polished.fastq|paste - - - -|sed 's/^@/>/'|awk '{print $1"\n"$2}' > pb_279_step3_v${i}_polished.fasta &  
  wait
  echo "Done" >> ../../countDone.txt ) & 
done 
for ((k=1;k<=$versions2;k++))
do
( cd ~/LUfungiProject/pb-assembly/pb_320-2_v${k}
  cat 2-asm-falcon/p_ctg.fa 2-asm-falcon/a_ctg.fa > pb_320-2_falcon_step2_v${k}.fasta
  cat 3-unzip/all_p_ctg.fa 3-unzip/all_h_ctg.fa >pb_320-2_falcon_step3_v${k}.fasta
  mv 2-asm-falcon/*.fa .
  mv 3-unzip/*.fa .
  find -maxdepth 1|grep \.\/[0-9].|while read dir;do rm -rf $dir;done
  ls ../../bamfiles/pb_320-2/*subreads.bam|while read bam;do echo $bam >> pb_320-2_bam.fofn;done
  nohup pbmm2 align pb_320-2_falcon_step3_v${k}.fasta pb_320-2_bam.fofn pb_320-2_v${k}_aligned.bam --sort -j 8 -J 8 -m 32G --preset SUBREAD & 
  nohup samtools faidx pb_320-2_falcon_step3_v${k}.fasta -o pb_320-2_falcon_step3_v${k}.fasta.fai &
  wait
  nohup pbindex pb_320-2_v${k}_aligned.bam &
  wait
  nohup arrow pb_320-2_v${k}_aligned.bam -r pb_320-2_falcon_step3_v${k}.fasta -o pb_320-2_step3_v${k}_polished.fastq -j 20 --diploid &
  wait
  nohup cat pb_320-2_step3_v${k}_polished.fastq|paste - - - -|sed 's/^@/>/'|awk '{print $1"\n"$2}' > pb_320-2_step3_v${k}_polished.fasta & 
  wait 
  echo "Done" >> ../../countDone.txt ) &
done 
cd ~/LUfungiProject/
sum=$(($versions1+$versions2))
count=$(cat countDone.txt|wc -l)
while [ ! "$count" == "$sum" ]
do 
  sleep 10
  count=$(cat countDone.txt|wc -l)
done
rm ~/LUfungiProject/countDone.txt

# Install busco and do busco analysis instead
conda create -n your_env_name -c bioconda -c conda-forge busco=4.0.6 python=3.7
conda activate busco
cd ~/LUfungiProject/pb-assembly/pb_279_v1
nohup busco -m genome -i pb_279_falcon_step3_v1.fasta -o step3_busco -l fungi_odb10 & #https://busco.ezlab.org/frames/fungi.html
cd ~/LUfungiProject/pb-assembly/pb_279_v2
nohup busco -m genome -i pb_279_falcon_step3_v2.fasta -o step3_busco -l fungi_odb10 &
cd ~/LUfungiProject/pb-assembly/pb_320-2_v1
nohup busco -m genome -i pb_320-2_falcon_step3_v1.fasta -o step3_busco -l fungi_odb10 &
cd ~/LUfungiProject/pb-assembly/pb_320-2_v2
nohup busco -m genome -i pb_320-2_falcon_step3_v2.fasta -o step3_busco -l fungi_odb10 &
# Do summary plotting for each busco analysis
mkdir /home2/shared_bioinformatics_master_projects/agaricalesGenomes/jiawei_zhao_assemblies/busco_plots
cd ~/LUfungiProject/pb-assembly/pb_?_v?/step3_busco
mkdir summaries
cp short_summary.specific.fungi_odb10.step3_busco.txt summaries/
nohup generate_plot.py -wd summaries/ &
cp summaries/*.png /home2/shared_bioinformatics_master_projects/agaricalesGenomes/jiawei_zhao_assemblies/busco_plots

#create a more automated analysis workflow for busco analysis 
touch ~/LUfungiProject/countDone.txt
cd ~/LUfungiProject/pb-assembly/
conda create -n busco -c bioconda -c conda-forge busco=4.0.6 python=3.7
conda deactivate
conda activate busco
versions1=$(find -maxdepth 1 -name "pb_279_v*"|wc -l);
versions2=$(find -maxdepth 1 -name "pb_320-2_v*"|wc -l);
for ((i=1;i<=$versions1;i++))
do
( cd ~/LUfungiProject/pb-assembly/pb_279_v${i}
  nohup busco -m genome -i pb_279_falcon_step2_v${i}.fasta -o step2_busco -l fungi_odb10 & sleep 10
  nohup busco -m genome -i pb_279_falcon_step3_v${i}.fasta -o step3_busco -l fungi_odb10 & sleep 10
  nohup busco -m genome -i pb_279_step3_v${i}_polished.fasta -o step4_busco -l fungi_odb10 & sleep 10
  wait
  echo "Done" >> ../../countDone.txt ) &
done 
for ((k=1;k<=$versions2;k++))
do
( cd ~/LUfungiProject/pb-assembly/pb_320-2_v${k}
  nohup busco -m genome -i pb_320-2_falcon_step2_v${k}.fasta -o step2_busco -l fungi_odb10 & sleep 10
  nohup busco -m genome -i pb_320-2_falcon_step3_v${k}.fasta -o step3_busco -l fungi_odb10 & sleep 10
  nohup busco -m genome -i pb_320-2_step3_v${k}_polished.fasta -o step4_busco -l fungi_odb10 & sleep 10
  wait
  echo "Done" >> ../../countDone.txt ) & 
done 
cd ~/LUfungiProject/
sum=$(($versions1+$versions2))
count=$(cat countDone.txt|wc -l)
while [ ! "$count" == "$sum" ]
do 
  sleep 10
  count=$(cat countDone.txt|wc -l)
done
rm ~/LUfungiProject/countDone.txt

# Generalize the busco-plotting pipeline
cd ~/LUfungiProject/pb-assembly/
mkdir busco_plots
cd busco_plots
ls ../*/step*_busco/*.txt|while read txt;do name=$(echo $txt|cut -d / -f 2-3|tr -d \/|sed 's/step/_step/');echo $name;mkdir $name;cp $txt $name;done
nohup sh -c 'ls *busco/|grep ^pb|sed 's/://'|while read dir;do nohup generate_plot.py -wd $dir;done' & wait
rm *.log
mv *_busco/ /home2/shared_bioinformatics_master_projects/agaricalesGenomes/jiawei_zhao_assemblies/busco_plots/
# Run busco analysis for the reference genomes
cd ~/LUfungiProject/OriginalAssemblies/
nohup busco -m genome -i pb_279_Leuge.fasta -o pb_279_ref_busco -l fungi_odb10 &
nohup busco -m genome -i pb_320-2_Mysco.fasta -o pb_320-2_ref_busco -l fungi_odb10 &
wait
find -maxdepth 2 -name "*busco.txt"|while read txt;do dir=$(echo $txt|cut -d / -f 1-2);mkdir ${dir}/summary;mv $txt ${dir}/summary;nohup generate_plot.py -wd ${dir}/summary/;done 
find -name "*summary"|while read dir;do name=$(echo $dir|cut -d / -f 2);pardir=$(echo $dir|cut -d / -f 1-2);newname=$(echo ${pardir}/$name);mv $dir $newname;done
mv */*_busco/ /home2/shared_bioinformatics_master_projects/agaricalesGenomes/jiawei_zhao_assemblies/busco_plots/

# Do quast analysis for assemblies
cd ~/LUfungiProject/pb-assembly/pb_279_v2
nohup quast.py -o step2_v2_quast/  pb_279_falcon_step2_v2.fasta -r ../../OriginalAssemblies/pb_279_Leuge.fasta -t 20 &
nohup quast.py -o step3_v2_quast/  pb_279_falcon_step3_v2.fasta -r ../../OriginalAssemblies/pb_279_Leuge.fasta -t 20 &
cd ~/LUfungiProject/pb-assembly/pb_320-2_v2
nohup quast.py -o step2_v2_quast/ pb_320-2_falcon_step2_v2.fasta -r ../../OriginalAssemblies/pb_320-2_Mysco.fasta -t 20 &
nohup quast.py -o step3_v2_quast/ pb_320-2_falcon_step3_v2.fasta -r ../../OriginalAssemblies/pb_320-2_Mysco.fasta -t 20 &

#Create a more automated analysis workflow for quast analysis
touch ~/LUfungiProject/countDone.txt
cd ~/LUfungiProject/pb-assembly/
conda deactivate
conda activate py2
versions1=$(find -maxdepth 1 -name "pb_279_v*"|wc -l);
versions2=$(find -maxdepth 1 -name "pb_320-2_v*"|wc -l);
for ((i=1;i<=$versions1;i++))
do
( cd ~/LUfungiProject/pb-assembly/pb_279_v${i}
  nohup quast.py -o step2_v${i}_quast/  pb_279_falcon_step2_v${i}.fasta -r ../../OriginalAssemblies/pb_279_Leuge.fasta -t 20 &
  nohup quast.py -o step3_v${i}_quast/  pb_279_falcon_step3_v${i}.fasta -r ../../OriginalAssemblies/pb_279_Leuge.fasta -t 20 &
  nohup quast.py -o step4_v${i}_quast/  pb_279_step3_v${i}_polished.fasta -r ../../OriginalAssemblies/pb_279_Leuge.fasta -t 20 & 
  wait
  echo "Done" >> ../../countDone.txt ) &
done 
for ((k=1;k<=$versions2;k++))
do
( cd ~/LUfungiProject/pb-assembly/pb_320-2_v${k}
  nohup quast.py -o step2_v${k}_quast/ pb_320-2_falcon_step2_v${k}.fasta -r ../../OriginalAssemblies/pb_320-2_Mysco.fasta -t 20 &
  nohup quast.py -o step3_v${k}_quast/ pb_320-2_falcon_step3_v${k}.fasta -r ../../OriginalAssemblies/pb_320-2_Mysco.fasta -t 20 &
  nohup quast.py -o step4_v${k}_quast/ pb_320-2_step3_v${k}_polished.fasta -r ../../OriginalAssemblies/pb_320-2_Mysco.fasta -t 20 & 
  wait
  echo "Done" >> ../../countDone.txt ) &
done 
cd ~/LUfungiProject/
sum=$(($versions1+$versions2))
count=$(cat countDone.txt|wc -l)
while [ ! "$count" == "$sum" ]
do 
  sleep 10
  count=$(cat countDone.txt|wc -l)
done
rm ~/LUfungiProject/countDone.txt

# Send my falcon assemblies and quast reports to the shared folder
cd ~/LUfungiProject/pb-assembly/
find -maxdepth 3 -name "report.pdf"|while read pdf;do dir=$(echo $pdf|cut -d / -f 1-3|sed -r 's/quast/quast\//');newname=$(echo $dir|cut -d / -f 2-3|sed -r 's/v[0-9]_//'|tr / _);mv "$pdf" "$dir${newname}.pdf";done
find -maxdepth 3 -name "*quast.pdf" #to check if pdf names are changed properly
#copy my quast reports and .fasta assemblies to the specific subfolder of the shared folder
find -maxdepth 3 -name "*quast.pdf"|while read pdf;do cp $pdf ../../../shared_bioinformatics_master_projects/agaricalesGenomes/jiawei_zhao_assemblies/pb-assembly/;done
find -maxdepth 2 -name "pb*.fasta"|while read asm;do cp $asm ../../../shared_bioinformatics_master_projects/agaricalesGenomes/jiawei_zhao_assemblies/pb-assembly/;done

# Run busco analysis for my previous assemblies
cd /home2/shared_bioinformatics_master_projects/agaricalesGenomes/jiawei_zhao_assemblies
find -mindepth 3 -name "*.gfa"|while read gfa;do name=$(echo $gfa|sed 's/.gfa//'|cut -d / -f 4);path=$(echo $gfa|cut -d / -f 1-3);nohup awk -v gfa="$gfa" -v name="$name" -v path="$path"  '/^S/{print ">"$2"\n"$3}' $gfa > ${path}/${name}.fasta;done #convert .gfa files to .fasta files
find -mindepth 3 -name "*.fasta"|while read fa;do (dir=$(echo $fa|cut -d / -f 1-3);fa=$(echo $fa|cut -d / -f 4);name=$(echo $fa|cut -d . -f 1);cd $dir;nohup busco -m genome -i $fa -o ${name}_busco -l fungi_odb10) & done
#However, busco analysis for some of those .fasta files failed since their headers consist of "/", which are forbidden for the input files of busco
#clean the headers of these .fasta files
ls */*/*.fasta|while read fasta;do fa=$(echo $fasta|sed 's/fasta/fa/');cat $fasta|tr -d "/" > $fa;done  
find -mindepth 3 -name "*.fa"|while read fa;do (dir=$(echo $fa|cut -d / -f 1-3);fa=$(echo $fa|cut -d / -f 4);name=$(echo $fa|cut -d . -f 1);cd $dir;nohup busco -m genome -i $fa -o ${name}_busco -l fungi_odb10) & done
#This command isn't gonna work when there're multiple busco tasks within one folder to be run simultaneously, manual operations should be done in order to prevent empty output folders
find -maxdepth 4 -mindepth 4 -name "*busco.txt"|while read txt;do dir=$(echo $txt|cut -d / -f 1-4);name=$(echo $txt|cut -d / -f 5|cut -d . -f 4);mkdir ${dir}/${name}_summary;mv $txt ${dir}/${name}_summary;nohup generate_plot.py -wd ${dir}/${name}_summary/;done
rm *.log;rm nohup.out
find -maxdepth 4 -mindepth 4 -name "*summary"|while read summary;do mv $summary ./busco_plots/;done
find -maxdepth 2 -name "*filtered*"|while read raven;do newname=$(echo $raven|sed 's/filtered_subreads/raven_/');mv $raven $newname;done
find -maxdepth 2 -name "*step*"|while read falcon;do newname=$(echo $falcon|sed 's/step/falcon_step/');mv $falcon $newname;done
find -maxdepth 2 -name "*falcon_falcon*"|while read falcon;do newname=$(echo $falcon|sed 's/_falcon_falcon_/_falcon_/');mv $falcon $newname;done
find -maxdepth 2 -name "*polished*summary"|while read miniasm;do newname=$(echo $miniasm|sed 's/polished_/polished_miniasm_/');mv $miniasm $newname;done
# In order to make the names of folders more recognizable

mv busco_plots/pb_279_busco_summary/ busco_plots/pb_279_flye_busco_summary/
mv busco_plots/pb_320-2_busco_summary/ busco_plots/pb_320-2_flye_busco_summary/



# Create a parallel version of quast analysis of the non-falcon assemblies
cd /home2/shared_bioinformatics_master_projects/agaricalesGenomes/jiawei_zhao_assemblies
mkdir quast_279
mkdir quast_320-2
find -mindepth 3 -name "*279*.fasta"|while read fa;do (dir=$(echo $fa|cut -d / -f 1-3);fasta=$(echo $fa|cut -d / -f 4);name=$(echo $fasta|cut -d . -f 1);nohup quast.py -o quast_279/${name}/ $fa -r ~/LUfungiProject/OriginalAssemblies/pb_279_Leuge.fasta -t 20)& done 
find -mindepth 3 -name "*320-2*.fasta"|while read fa;do (dir=$(echo $fa|cut -d / -f 1-3);fasta=$(echo $fa|cut -d / -f 4);name=$(echo $fasta|cut -d . -f 1);nohup quast.py -o quast_320-2/${name}/ $fa -r ~/LUfungiProject/OriginalAssemblies/pb_320-2_Mysco.fasta -t 20)& done #compare each of my assembly to the reference assembly
ls */*/report.pdf|while read report;do newname=$(echo $report|cut -d / -f 2);dir=$(echo $report|cut -d / -f 1);mv $report $dir/${newname}.pdf;done #change the names of the .pdf reports in order to make them recognizable

#Generate fastqc reports for the output files of .fastq format from genomicconsensus
cd ~/LUfungiProject/pb-assembly/
mkdir fastqcReports;ls */*polished.fastq|while read fastq;do (name=$(echo $fastq|cut -d / -f 2|cut -d . -f 1);mkdir fastqcReports/$name;nohup fastqc --threads 20 -o fastqcReports/$name -f fastq $fastq ) & done
#Unfortunately, fastqc analysis of two .fastq files failed with "java:out of memory exception"

#multiqc1.9 becomes installed and available to analyse busco & quast results
cd ~/LUfungiProject/pb-assembly/
conda deactivate
conda activate py3.6
multiqc busco_plots279/ -o multiqc_busco279 &>> mtqc.log &
multiqc busco_plots320-2/ -o multiqc_busco320-2 &>> mtqc.log &
multiqc pb_320-2_*/*quast/report.tsv -o multiqc_quast320-2 &>> mtqc.log &
multiqc pb_279_*/*quast/report.tsv -o multiqc_quast279 &>> mtqc.log &


