#This script records bash commands used for fungal genome assembly

##acquire information about all subreads and filtered subreads in the dataset
ls -lah ../../shared_bioinformatics_master_projects/agaricalesGenomes/b2016040/INBOX/*/*subreads/;
##create correspnding repositories in the rawdata  directory
ls ../../shared_bioinformatics_master_projects/agaricalesGenomes/b2016040/INBOX/*/*subreads/*|cut -d / -f 7|while read fungi_name;do mkdir ./rawdata/$fungi_name;done;
##Do fastqc analysis for the 4 given fastq.gz files
ls fastqcReports/*|cut -d / -f 2|sed s/://g|while read name;do nohup fastqc --threads 20 -o ./fastqcReports// --extract -f fastq -c ../../shared_bioinformatics_master_projects/agaricalesGenomes/b2016040/INBOX//*reads/*.gz;done;
