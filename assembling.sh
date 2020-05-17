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
nohup ./wtdbg2 -x rs -g 6.5g -i ../pb_279_filtered_subreads.fastq.gz -t 20 -fo ../wtdbg2Assembly/pb_279_default/pb_279_default;
nohup ./wtdbg2 -x rs -g 6.9g -i ../pb_320-2_filtered_subreads.fastq.gz -t 20 -fo ../wtdbg2Assembly/pb_320-2_default/pb_320-2_default;
# derive raw consensus
ls ../wtdbg2Assembly/*/*.ctg.lay.gz|while read path;do dir=$(echo $path|cut -d / -f 1-3);name=$(echo $path|cut -d / -f 4|cut -d . -f 1);nohup ./wtpoa-cns -t 20 -i $path -fo $dir/$name.raw.fa;done;

##In order to convert original .fastq.gz files of short reads into .BAM files for genome polishing with arrow, it's necessary to use blasr(bowtie2 is used for short reads, inappropriate here)
ls ../unzipped/*.fastq|while read subreads;do name=$(echo $subreads|cut -d / -f 3|cut -d _ -f 1,2);ref=$(echo ./$name*/$name*.raw.fa);dir=$(echo $ref|cut -d / -f 1,2);nohup blasr $subreads $ref --sam --out $dir/blasrBAM/$name.sam;done  #however, blasr doesn't work here, 'your subreads aren't PacBio reads'

##Long read assembly by Miniasm
# Find overlaps by all-vs-all self-mappings via minimap2
nohup minimap2 -x ava-pb -t 15 pb_279_filtered_subreads.fastq.gz pb_279_filtered_subreads.fastq.gz|gzip -1 >miniasmAssembly/pb_279_default.paf.gz
nohup minimap2 -x ava-pb -t 15 pb_320-2_filtered_subreads.fastq.gz pb_320-2_filtered_subreads.fastq.gz|gzip -1 >miniasmAssembly/pb_320-2_default.paf.gz
# Layout using Miniasm
ls miniasmAssembly/*|while read paf;do name=$(echo $paf|cut -d / -f 2|cut -d . -f 1|cut -d _ -f 1,2);nohup miniasm -f $name*.fastq.gz $paf > ./miniasmAssembly/$name.gfa;done 
# Convert .gfa to .fasta
awk '/^S/{print ">"$2"\n"$3}' pb_279.gfa|fold > pb_279-unpolished.fasta
awk '/^S/{print ">"$2"\n"$3}' pb_320-2.gfa|fold > pb_320-2_unpolished.fasta

# Assembly polishing by minipolish(the output files aren't actually in .gfa format, eerie here)
nohup minipolish -t 20 --pacbio pb_279_filtered_subreads.fastq.gz miniasmAssembly/pb_279.gfa > miniasmAssembly/pb_279_polished.gfa
nohup minipolish -t 20 --pacbio pb_320-2_filtered_subreads.fastq.gz miniasmAssembly/pb_320-2.gfa > miniasmAssembly/pb_320-2_polished.gfa

## Genome assembling and polishing by raven with default setting
ls *.fastq.gz|while read fastq;do name=$(echo $fastq|cut -d . -f 1);nohup raven --graphical-fragment-assembly ravenAssembly/default/$name.gfa  -t 20  $fastq;echo $name;done 
# Convert .gfa files to .fasta files
awk '/^S/{print ">"$2"\n"$3}' pb_279_filtered_subreads.gfa|fold >pb_279.fasta;
awk '/^S/{print ">"$2"\n"$3}' pb_320-2_filtered_subreads.gfa|fold >pb_320-2.fasta; 
