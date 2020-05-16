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
# derive consensus
ls ../wtdbg2Assembly/*/*.ctg.lay.gz|while read path;do dir=$(echo $path|cut -d / -f 1-3);name=$(echo $path|cut -d / -f 4|cut -d . -f 1);nohup ./wtpoa-cns -t 20 -i $path -fo $dir/$name.raw.fa;done;

##In order to convert original .fastq.gz files of short reads into .BAM files for genome polishing with arrow, it's necessary to use blasr(bowtie2 is used for short reads, inappropriate here)
ls ../unzipped/*.fastq|while read subreads;do name=$(echo $subreads|cut -d / -f 3|cut -d _ -f 1,2);ref=$(echo ./$name*/$name*.raw.fa);dir=$(echo $ref|cut -d / -f 1,2);nohup blasr $subreads $ref --sam --out $dir/blasrBAM/$name.sam;done  #however, there were errors within this step and the output files seem like incomplete

##Try to use pbalign instead

