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
ls -d */|sed "s/\///g"|while read name;do echo $name >> .gitignore;done

##Long PacBio RSII reads assembly by wtdbg2 
nohup ./wtdbg2 -x rs -g 6.5g -i ../pb_279_filtered_subreads.fastq.gz -t 20 -fo ../wtdbg2Assembly/pb_279_default/pb_279_default
nohup ./wtdbg2 -x rs -g 6.9g -i ../pb_320-2_filtered_subreads.fastq.gz -t 20 -fo ../wtdbg2Assembly/pb_320-2_default/pb_320-2_default

##In order to convert original .fastq.gz files of short reads into .BAM files for genome polishing with arrow, it's necessary to first use bowtie2 to have SAM format which can be transformed to BAM format with Samtools then
 
