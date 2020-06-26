<h1 align="center">TWO GENOME ASSEMBLY PROJECT</h1>
<p align="left">Litter decomposing fungi represent important agents of soil carbon cycling across terrestrial ecosystems ranging from tundra to tropical forests. The reasons behind the enormous success of litter decomposers are not well known, but it is reasonable to consider that this could be partly due to the diverse plant cell wall decomposition systems these fungi harbor. To examine the diversity of these systems in litter decomposers, we sequenced nine species with diverse habitat preferences from Agaricales, which represents one of the most species-rich order in mushroom forming fungi and in which many litter decomposers are found. All genomes were sequenced using PacBio and the assemblies of seven of these genomes were of good quality. However, the genome assemblies of Mycetinis scorodonius(Mysco, pb_320-2) and Leucopaxillus gentianeus(Leuge, pb_279) did not pass the strict quality standards, that we have set up in this project. The reasons for this could be related to the larger than expected genome size for both species in combination to their dikaryotic state. The successful genome assembly for the two species is of particular interest since both species are found exclusively on conifer litter and therefore, may harbor adaptations related to this recalcitrant type of substrate. Genome assembly tools for dikaryotic genomes are continuously improving, which poses a great opportunity to attempt again the assembly of the two remaining genomes using recently developed tools. </p>

***
[**Rawdata-processing tools**](#intro0)
+ [**bax2bam**](#bax2bam)
***
[**Genome assemblers and polishers**](#intro1)
+ [**pb-assembly**](#pbasm)
+ [**pbmm2**](#pbmm2)
+ [**genomicconsensus**](#gcc)
+ [**miniasm**](#miniasm)
+ [**minipolish**](#minipolish)
+ [**raven**](#raven)
+ [**flye**](#flye)
+ [**canu**](#canu)
***
[**Evaluation tools**](#intro2)
+ [**fastqc**](#fastqc)
+ [**multiqc**](#multiqc)
+ [**QUAST**](#quast)
+ [**Bandage**](#bandage)
+ [**BUSCO**](#busco)
***
[**Pipelines of analysis**](#pipelines)
+ [**Preliminary data-analysis & data-processing**](#preli)
+ [**Try different non-falcon long-read assemblers to generate genome assemblies**](#trydi)
+ [**Use pb-assembly(falcon and falcon-unzip) & genomicconsensus for assembling & polishing**](#falgen)
+ [**Assembly evaluation by QUAST & BUSCO**](#quabus)
***

<a name="intro0"></a>
# Introduction to software which would be used for preliminary data processing

<a name="bax2bam"></a>
## bax2bam
bax2bam converts the legacy PacBio basecall format (bax.h5) into the BAM basecall format.

### Installation and help
```bash
conda create -n py3 python=3.7
conda activate py3
conda install -c bioconda bax2bam
```
### Official Github Page
- [bax2bam](https://github.com/pacificbiosciences/bax2bam/)

<a name="intro1"></a>
# Introduction to software which would be used or tried for assembling genomes
Because of the dikaryotic nature of these two fungi, a chronological combination of falcon & falcon-unzip & falcon-phase from  [***pb-assembly***](#pbasm) is supposed to be the most appropriate combo here. Nonetheless, several other popular genome assemblers & polishers ([***miniasm***](#miniasm),[***minipolish***](#minipolish),[***raven***](#raven),[***flye***](#flye),[***canu***](#canu)) which are available for PacBio subreads would be listed and tried as well. (However, [***CANU***](#canu) assembler was ultimately discarded since it had already taken too much time and still hadn't given final assemblies).

Since I have encountered the same issue when using falcon-unzip as the unsolved issue here shows, 

<a name="error"></a>
## An unsolved error
- [Error at 4-polish](https://github.com/PacificBiosciences/FALCON_unzip/issues/159)

<a name="solution"></a>
### A compromised solution
I decided to concatenate the .fasta file of primary contigs and the .fasta file of haplotype contigs into a reference assembly file. Both the primary-contig file and the haplotype-contig file were output from the successful '3-unzip' step of falcon-unzip. Thus, the concatenated .fasta file would be used as a reference input for [***pbmm2***](#pbmm2) & for + [***genomicconsensus***](#gcc) in order to generate a polished consensus .fasta file. Unfortunately, it would be impossible to use falcon-phase then. As a result, the consensus .fasta file would be considered as the final assembly. In order to access more detailed descriptions of PacBio SMRT tools, an online tutorial is available at [***SMRTtools***](https://www.pacb.com/wp-content/uploads/SMRT_Tools_Reference_Guide_v700.pdf).  

<a name="pbasm"></a>
## pb-assembly
pb-assembly is an official PacBio combo of three genome assembly tools - falcon, falcon-unzip and falcon-phase. It provides an integrated workflow to assemble diploid genomes from subread files whose formats are .fasta and .bam.
### Installation and help

```bash
conda create -n py2 python=2.7
conda activate py2
conda install -c bioconda/label/cf201901 nim-falcon #https://github.com/PacificBiosciences/pbbioconda/issues/265
conda install -c bioconda pb-assembly
conda install -c bioconda samtools=1.9 --force-reinstall #https://github.com/bioconda/bioconda-recipes/issues/13958
fc_run.py -h
```
### Official Github Page
- [PacBio:falcon](https://github.com/PacificBiosciences/pb-assembly)

<a name="pbmm2"></a>
## pbmm2
pbmm2 is a SMRT C++ wrapper for minimap2's C API. Its purpose is to support native PacBio in- and output, provide sets of recommended parameters, generate sorted output on-the-fly, and postprocess alignments. Sorted output can be used directly for polishing using GenomicConsensus, if BAM has been used as input to pbmm2. Benchmarks show that pbmm2 outperforms BLASR in sequence identity, number of mapped bases, and especially runtime. pbmm2 is the official replacement for BLASR. It would require an .fofn file which consists of paths of unaligned .bam subread files and a reference genome as input.     

### Installation and help
```bash
conda activate py2
conda install -c bioconda pbmm2
pbmm2 -h
```
### Official Github Page
- [PacBio:pbmm2](https://github.com/PacificBiosciences/pbmm2)

<a name="gcc"></a>
## GenomicConsensus
GenomicConsensus package provides a variant-calling tool which allows users to apply one of the algorithms within {quiver,arrow,plurality,poa,best} to compute genomic consensus and call variants relative to the reference. An aligned .bam file and a reference genome are required.

### Installation and help

```bash
conda activate py2
conda install -c bioconda genomicconsensus python=2
quiver -h
arrow -h
```
### Official Github Page
- [PacBio:genomicconsensus](https://github.com/PacificBiosciences/GenomicConsensus)

<a name="miniasm"></a>
## miniasm
Miniasm is a rather particular long-read assembler as it doesn't include a consensus step. The resulting contigs are just merged erroneous long reads and still contains many sequencing errors. Produced contigs are structually correct, but at the nucleotide level, there're many mismatches and indels. Its output format is .gfa.

### Installation and help

```bash
conda activate py2
conda install -c bioconda miniasm
miniasm -h
minimap2 -h
```

### Official Github Page
- [miniasm](https://github.com/lh3/miniasm)

<a name="minipolish"></a>
## minipolish
Minipolish is a combined tool for .gfa genome assembly polishing, especially for those which were produced by miniasm. It combines racon and minimap and can run for multiple rounds. It supports both ONT and PACBIO subreads as input.

### Installation and help

```bash
conda activate py3
conda install -c bioconda minipolish
minipolish -h
```
### Official github page
- [minipolish](https://github.com/rrwick/Minipolish#installation)

<a name="raven"></a>
## raven
Raven is an advanced version of Racon-assembler. It's a de novo genome assembler for long uncorrected reads. It's likely to be relatively reliable for chromosome assembly, though it doesn't perform well on small plasmids and has circularisation issues.

### Installation and help

```bash
conda activate py3
conda install -c bioconda raven-assembler 
raven -h
```
### Official github page
- [raven-assembler](https://github.com/lbcb-sci/raven)

<a name="flye"></a>
## flye
Flye is a de novo assembler for single molecule sequencing reads, such as those produced by PacBio and Oxford Nanopore Technologies. It is designed for a wide range of datasets, from small bacterial projects to large mammalian-scale assemblies. The package represents a complete pipeline: it takes raw PacBio / ONT reads as input and outputs polished contigs. Flye also has a special mode for metagenome assembly.

### Installation and help

```bash
conda activate py3
conda install -c bioconda flye
flye -h
```
### Official github page
- [flye-assembler](https://github.com/fenderglass/Flye)

<a name="canu"></a>
## canu
Canu is a fork of the Celera Assembler, designed for high-noise single-molecule sequencing (such as the PacBio RS II/Sequel or Oxford Nanopore MinION).
Canu is a hierarchical assembly pipeline which runs in four steps:

- Detect overlaps in high-noise sequences using MHAP

- Generate corrected sequence consensus

- Trim corrected sequences

- Assemble trimmed corrected sequences

### Installation and help

```bash
conda activate py3
conda install -c bioconda canu
export PERL5LIB="" #https://groups.google.com/a/continuum.io/g/anaconda/c/5DFTW1GDgXQ?pli=1
canu -h
```
### Official github page
- [canu-assembler](https://github.com/marbl/canu)

<a name="intro2"></a>
# Introduction to software which would be used or tried for evaluation of raw data & assembly evaluation

<a name="fastqc"></a>
## fastqc
FastQC is a program designed to spot potential problems in high througput sequencing datasets. It runs a set of analyses on one or more raw sequence files in fastq or bam format and produces a report which summarises the results.

FastQC will highlight any areas where this library looks unusual and where you should take a closer look. The program is not tied to any specific type of sequencing technique and can be used to look at libraries coming from a large number of different experiment types (Genomic Sequencing, ChIP-Seq, RNA-Seq, BS-Seq etc etc).

### Installation and help

```bash
conda activate py2
conda install -c bioconda fastqc
fastqc -h
```
### Official github page
- [fastqc](https://github.com/s-andrews/FastQC)

<a name="multiqc"></a>
## multiqc
MultiQC is a tool to create a single report with interactive plots for multiple bioinformatics analyses across many samples.
There a very large number of Bioinformatics tools supported by MultiQC. Please see the MultiQC website for a [***complete list***](https://multiqc.info/#supported-tools).

MultiQC is written in Python (tested with v3.6+). It is available on the Python Package Index and through conda using Bioconda.

Reports are generated by scanning given directories for recognised log files. These are parsed and a single HTML report is generated summarising the statistics for all logs found. MultiQC reports can describe multiple analysis steps and large numbers of samples within a single plot, and multiple analysis tools making it ideal for routine fast quality control.

### Installation and help

```bash
conda create -n py3.6 -c conda-forge -c bioconda multiqc=1.9 python=3.6
conda activate py3.6
multiqc -h
```
### Official github page and official documentation
- [multiqc-github](https://github.com/ewels/MultiQC)
- [multiqc-documentation](https://multiqc.info/docs/#running-multiqc)

<a name="quast"></a>
## QUAST
QUAST stands for QUality ASsessment Tool. It evaluates genome/metagenome assemblies by computing various metrics. The current QUAST toolkit includes the general QUAST tool for genome assemblies, MetaQUAST, the extension for metagenomic datasets, QUAST-LG, the extension for large genomes (e.g., mammalians), and Icarus, the interactive visualizer for these tools.

The QUAST package works both with and without reference genomes. However, it is much more informative if at least a close reference genome is provided along with the assemblies. The tool accepts multiple assemblies, thus is suitable for comparison.

### Installation and help

```bash
conda activate py2
conda install -c bioconda quast
quast -h
```
### Official github page
- [QUAST](https://github.com/ablab/quast)

<a name="bandage"></a>
## Bandage
Bandage is a GUI program that allows users to interact with the assembly graphs made by de novo assemblers such as [***Velvet***](https://www.ebi.ac.uk/~zerbino/velvet/), [***SPAdes***](http://bioinf.spbau.ru/spades), [***MEGAHIT***](https://github.com/voutcn/megahit) and others.

De novo assembly graphs contain not only assembled contigs but also the connections between those contigs, which were previously not easily accessible. Bandage visualises assembly graphs, with connections, using graph layout algorithms. Nodes in the drawn graph, which represent contigs, can be automatically labelled with their ID, length or depth. Users can interact with the graph by moving, labelling and colouring nodes. Sequence information can also be extracted directly from the graph viewer. By displaying connections between contigs, Bandage opens up new possibilities for analysing and improving de novo assemblies that are not possible by looking at contigs alone.

### Installation and help

```bash
conda activate py3
conda install -c bioconda bandage
Bandage -h
```
### Official github page
- [bandage](https://github.com/rrwick/Bandage)


<a name="busco"></a>
## busco
BUSCO - Benchmarking Universal Single-Copy Orthologs.

Based on evolutionarily-informed expectations of gene content of near-universal single-copy orthologs, BUSCO metric is complementary to technical metrics like N50. 

Besides evaluation of genomes, BUSCO can be used for evaluation of proteins & of transcriptomes as well.

### Installation and help

```bash
conda create -n busco -c bioconda -c conda-forge busco=4.0.6 python=3.7
conda activate busco
busco -h
```
### Official github page and official documentation
- [busco-github](https://github.com/openpaul/busco)
- [busco-documentation](https://busco.ezlab.org/)

<a name="pipelines"></a>
# Pipelines of analysis
First of all, the file-tree of my local repository on the server is shown here: [***localRepoTree***](https://github.com/Orthologues/LUfungiProject/blob/master/LUfungiProject.tree). 

What's more, the absolute path of the local directory on the server which stores all rawdata is "/home2/shared_bioinformatics_master_projects/agaricalesGenomes". The file-tree of this rawdata directory is shown here: [***rawdataDirTree***](https://github.com/Orthologues/LUfungiProject/blob/master/rawdata.tree).

In order to achieve achieve availability of my integrated pipeline scripts of [**assembling by falcon & quast analysis**](https://github.com/Orthologues/LUfungiProject/blob/master/pb-assembly/integrated_falcon_quast.sh) and of [**falcon-assembly busco analysis**](https://github.com/Orthologues/LUfungiProject/blob/master/pb-assembly/integrated_busco.sh), your local repository of a reproduced diploid genome-assembly project must be fulfill the following requirements:

- Is directly under your own home("\~/") directory on the server(system). In my own case here, my local "LUfungiProject" repository has a relative path as "\~/LUfungiProject".

- Follows the exact subdirectory & file structure and naming formats of my local repository. 

<a name="preli"></a>
## Preliminary data-analysis & data-processing

```bash
cd ~/LUfungiProject

## FASTQC 
mkdir fastqcReports #Create a directory which stores fastqcReports of subread files
#create a corresponding repository for each of the 4 given fungal species in "fastqcReports" directory
ls ../../shared_bioinformatics_master_projects/agaricalesGenomes/b2016040/INBOX/*/*subreads/*|cut -d / -f 7|while read fungi_name;do mkdir ./fastqcReports/$fungi_name;done; 
#Do fastqc analysis for the given raw subread file in fastq.gz format of each of the 4 fungal species 
conda activate py2
nohup sh -c 'ls fastqcReports/*|cut -d / -f 2|sed s/://g|while read name;do nohup fastqc --threads 20 -o ./fastqcReports/$name/ --extract -f fastq -c ../../shared_bioinformatics_master_projects/agaricalesGenomes/b2016040/INBOX/$name/*reads/*.gz;done' & wait 
#generate a integrated multiqc report from all the preceding fastqc reports of .fastq.gz subread files
conda activate py3.6 #multiqc must be run under this env
nohup multiqc ./fastqcReports/ -o multiQCreport/ & wait 

## SUBREAD FILES
#Create symbolic links for the 2 subread files in .fastq.gz which would be necessary later
ln -s ../../shared_bioinformatics_master_projects/agaricalesGenomes/b2016040/INBOX/pb_279/filtered_subreads/pb_279_filtered_subreads.fastq.gz .
ln -s ../../shared_bioinformatics_master_projects/agaricalesGenomes/b2016040/INBOX/pb_320-2/subreads/pb_320-2_filtered_subreads.fastq.gz .  
#Convert .fastq.gz files to .fasta files as required
nohup less pb_279_filtered_subreads.fastq.gz|paste - - - -|sed 's/^@/>/'|awk '{print $1"\n"$2}' > pb_279_filtered_subreads.fasta &
nohup less pb_320-2_filtered_subreads.fastq.gz|paste - - - -|sed 's/^@/>/'|awk '{print $1"\n"$2}' > pb_320-2_filtered_subreads.fasta &
wait
#Calculate the sequence length distribution of .fasta files
nohup sh -c 'ls *.fasta|while read file;do nohup awk -v file="$file"  '/^[^>]/ {sum++;for(i=2000;i<=20000;i+=1000){if(length($0)>=i){qualified[(i-2000)/1000]++}}} END{print sum >> "lenstat.txt";for(i=0;i<=18;i++){printf("%d - %.2f\n",i*1000+2000,qualified[i]/sum*100) >> "lenstat.txt"}}' $file;done' & wait

## REFERENCE ASSEMBLIES
#Concatenate the primary-contig.fasta file and the associate-contig.fasta file of each of the two fungi into a reference assembly
mkdir OriginalAssemblies
cat /home2/shared_bioinformatics_master_projects/agaricalesGenomes/genome_assemblies/pb_279_Leuge/not_polished/*.fa >
OriginalAssemblies/pb_279_Leuge.fasta
cat /home2/shared_bioinformatics_master_projects/agaricalesGenomes/genome_assemblies/pb_320-2_Mysco/not_polished/*.fa >
OriginalAssemblies/pb_320-2_Mysco.fasta

## BAX2BAM
#Create symbolic links for the directories which contain .bax.h5 files which would be necessary to generate .bam files
ln -s ../../shared_bioinformatics_master_projects/agaricalesGenomes/b2016040/INBOX/pb_279/rawdata/run1 pb_279_raw
ln -s ../../shared_bioinformatics_master_projects/agaricalesGenomes/b2016040/INBOX/pb_320-2/raw/run1 pb_320-2_raw 
mkdir bamfiles #Create a directory which stores bamfiles converted by bax2bam from bax.h5 files
#Create a general input list for all .bax.h5 files of each of the two fungal species  
ls pb_279_raw/*/Analysis_Results/*.bax.h5 >> pb_279_baxh5.txt
ls pb_320-2_raw/*/Analysis_Results/*.bax.h5 >> pb_320-2_baxh5.txt 
#Since each of the 2 species have 24 .bax.h5 files which are from 8 movies and each of the 8 movies have 3 .bax.h5 files, 2*8 bam files would be generated. Since bax2bam needs an file of .bax.h5 file names for producing every .bam file, 2 for-loops are used here to produce 2*8 files of file names     
for i in {3..24..3};do sed -n "$(echo "$i-2"|bc),${i}p" < pb_279_baxh5.txt >> bamfiles/pb_279_list/pb_279_bax_list$(echo "$i/3"|bc).txt;done
for i in {3..24..3};do sed -n "$(echo "$i-2"|bc),${i}p" < pb_320-2_baxh5.txt >> bamfiles/pb_320-2_list/pb_320-2_bax_list$(echo "$i/3"|bc).txt;done 
#Use bax2bam to generate .bam files which would be necessary for pb-assembly(Falcon)
conda activate py3 
nohup sh -c 'for i in {1..8..1};do nohup bax2bam -f bamfiles/pb_279_list/pb_279_bax_list${i}.txt -o bamfiles/pb_279/${i} --subread;done' &
nohup sh -c 'for i in {1..8..1};do nohup bax2bam -f bamfiles/pb_320-2_list/pb_320-2_bax_list${i}.txt -o bamfiles/pb_320-2/${i} --subread;done' &
wait 
```

<a name="trydi"></a>
## Try different non-Falcon long-read assemblers to generate genome assemblies
```bash
cd ~/LUfungiProject

## MINIASM
conda activate py2
mkdir miniasmAssembly
#find overlaps by all-vs-all self-mappings via minimap2
nohup minimap2 -x ava-pb -t 15 pb_279_filtered_subreads.fastq.gz pb_279_filtered_subreads.fastq.gz|gzip -1 >miniasmAssembly/pb_279_default.paf.gz &
nohup minimap2 -x ava-pb -t 15 pb_320-2_filtered_subreads.fastq.gz pb_320-2_filtered_subreads.fastq.gz|gzip -1 >miniasmAssembly/pb_320-2_default.paf.gz &
wait
#Layout using Miniasm
nohup sh -c 'ls miniasmAssembly/*|while read paf;do name=$(echo $paf|cut -d / -f 2|cut -d . -f 1|cut -d _ -f 1,2);nohup miniasm -f $name*.fastq.gz $paf > ./miniasmAssembly/${name}_miniasm_unpolished.gfa;done' & wait
#Convert .gfa files to .fasta files
awk '/^S/{print ">"$2"\n"$3}' miniasmAssembly/pb_279_miniasm_unpolished.gfa > miniasmAssembly/pb_279_miniasm_unpolished.fasta
awk '/^S/{print ">"$2"\n"$3}' miniasmAssembly/pb_320-2_miniasm_unpolished.gfa > miniasmAssembly/pb_320-2_miniasm_unpolished.fasta
#Assembly polishing by minipolish(nohup shouldn't be used here since it would redirect stdout to the .gfa file)
conda activate py3
minipolish -t 15 --pacbio pb_279_filtered_subreads.fastq.gz miniasmAssembly/pb_279_miniasm_unpolished.gfa > miniasmAssembly/pb_279_miniasm_polished.gfa &
minipolish -t 15 --pacbio pb_320-2_filtered_subreads.fastq.gz miniasmAssembly/pb_320-2_miniasm_unpolished.gfa > miniasmAssembly/pb_320-2_miniasm_polished.gfa &
wait
# Convert .gfa to .fasta
awk '/^S/{print ">"$2"\n"$3}' miniasmAssembly/pb_279_miniasm_polished.gfa > miniasmAssembly/pb_279_miniasm_polished.fasta
awk '/^S/{print ">"$2"\n"$3}' miniasmAssembly/pb_320-2_miniasm_polished.gfa > miniasmAssembly/pb_320-2_miniasm_polished.fasta

## RAVEN
mkdir ravenAssembly
conda activate py3
#Genome assembling and polishing by raven with default setting
nohup sh -c 'ls *.fastq.gz|while read fastq;do name=$(echo $fastq|cut -d . -f 1|cut -d _ -f 1-2);nohup raven --graphical-fragment-assembly ravenAssembly/${name}_raven_default.gfa  -t 20 $fastq;done' & wait 
#change the alignment parameters
#match score=2, mismatch penalty=-6, gap penalty=-4, polishing iterations=2
nohup sh -c 'ls *.fastq.gz|while read fastq;do name=$(echo $fastq|cut -d . -f 1|cut -d _ -f 1-2);nohup raven -m 2 -n -6 -g -4  --graphical-fragment-assembly ravenAssembly/${name}_raven_m2n-6g-4p2.gfa  -t 20 $fastq;done' & wait 
#match score=3, mismatch penalty=-5, gap penalty=-4, polishing iterations=3
nohup sh -c 'ls *.fastq.gz|while read fastq;do name=$(echo $fastq|cut -d . -f 1|cut -d _ -f 1-2);nohup raven -p 3 -m 3 -n -5 -g -4 --graphical-fragment-assembly ravenAssembly/${name}_raven_m3n-5g-4p3.gfa -t 20 $fastq;done' & wait 
#match score=2, mismatch penalty=-6, gap penalty=-4, polishing iterations=3
nohup sh -c 'ls *.fastq.gz|while read fastq;do name=$(echo $fastq|cut -d . -f 1|cut -d _ -f 1-2);nohup raven -p 3 -m 2 -n -6 -g -4 --graphical-fragment-assembly ravenAssembly/${name}_raven_m2n-6g-4p3.gfa -t 20 $fastq;done' & wait 
#match score=2, mismatch penalty=-6, gap penalty=-4, polishing iterations=1
nohup sh -c 'ls *.fastq.gz|while read fastq;do name=$(echo $fastq|cut -d . -f 1|cut -d _ -f 1-2);nohup raven -p 1 -m 2 -n -6 -g -4 --graphical-fragment-assembly ravenAssembly/${name}_raven_m2n-6g-4p1.gfa -t 20 $fastq;done' & wait 
#match score=3, mismatch penalty=-5, gap penalty=-4, polishing iterations=1
nohup sh -c 'ls *.fastq.gz|while read fastq;do name=$(echo $fastq|cut -d . -f 1|cut -d _ -f 1-2);nohup raven -p 1 -m 3 -n -5 -g -4 --graphical-fragment-assembly ravenAssembly/${name}_raven_m3n-5g-4p1.gfa -t 20 $fastq;done' & wait 
#convert all .gfa files into .fasta files
ls ravenAssembly/*.gfa|while read gfa;do name=$(echo $gfa|cut -d . -f 1);awk '/^S/{print ">"$2"\n"$3}' $gfa > ${name}.fasta;done

## FLYE
mkdir flyeAssembly
#'genome-size' here are set according to the sizes of corresponding reference assemblies 
nohup flye --pacbio-raw pb_279_filtered_subreads.fastq.gz --genome-size 74m --out-dir flyeAssembly/pb_279_default/ --threads 15 &
nohup flye --pacbio-raw pb_320-2_filtered_subreads.fastq.gz --genome-size 137m --out-dir flyeAssembly/pb_320-2_default/ --threads 15 &
wait
#rename and move flye assemblies
mv flyeAssembly/pb_279_default/assembly.fasta flyeAssembly/pb_279_flye_default.fasta
mv flyeAssembly/pb_279_default/assembly.gfa flyeAssembly/pb_279_flye_default.gfa
mv flyeAssembly/pb_320-2_default/assembly.fasta flyeAssembly/pb_320-2_flye_default.fasta
mv flyeAssembly/pb_320-2_default/assembly.gfa flyeAssembly/pb_320-2_flye_default.gfa

## CANU
#'genome-size' here are set according to the sizes of corresponding reference assemblies 
nohup canu -p pb_279 -d canuAssembly/pb_279_default/ genomeSize=74m errorRate=0.3 gnuplotTested=true -pacbio-raw pb_279_filtered_subreads.fastq.gz 
nohup canu -p pb_320-2 -d canuAssembly/pb_320-2_default/ genomeSize=137m errorRate=0.3 gnuplotTested=true -pacbio-raw pb_320-2_filtered_subreads.fastq.gz 
#Unfortunately, CANU had taken too much time and memories so I decided to abort it
pkill canu
rm -rf canuAssembly
```

<a name="falgen"></a>
## Use pb-assembly(falcon and falcon-unzip) & genomicconsensus for assembling & polishing
```bash
cd ~/LUfungiProject
conda activate py2
mkdir pb-assembly
#Create a repository which stores configuration files of falcon-assembler(fc_run.py) and falcon-unzip(fc_unzip.py) & falcon-phase(fc_phase.py)
mkdir pb-assembly/mycfgs/
#Create .fofn files which contain paths of .fasta and .bam files for falcon & falcon-unzip
cd pb-assembly
ls ../bamfiles/pb_279/*subreads.bam|while read bam;do echo $bam >> pb_279_bam.fofn;done
ls ../bamfiles/pb_320-2/*subreads.bam|while read bam;do echo $bam >> pb_320-2_bam.fofn;done
echo "../pb_279_filtered_subreads.fasta" >> pb_279_fa.fofn 
echo "../pb_320-2_filtered_subreads.fasta" >> pb_320-2_fa.fofn 
```
Paths of the .fofn files which were created above would be filled into my configuration files as one of their compulsory parameters. My configuration files can be viewed [***here***](https://github.com/Orthologues/LUfungiProject/tree/master/pb-assembly/mycfgs). 'v1' means version1, 'v2' means version2, etc.
```bash
#Create a corresponding working directoriy for each of the fc_pb* configuration files.
for i in {1..2};do mkdir pb_279_v${i}/;done
for i in {1..3};do mkdir pb_320-2_v${i}/;done
```
There are multiple working directories and the pipeline of assembly from each .cfg file is relatively complicated. Therefore, I would like to use only one of my working directories - 'pb_279_v2' as an example to show the pipeline of assembly
```bash
cd pb_279_v2
#run falcon-assembly
nohup fc_run ../mycfgs/fc_pb_279_v2.cfg &> run0.log & wait 
#run falcon-unzip
nohup fc_unzip.py ../mycfgs/fc_unzip_pb_279.cfg &> run1.std &
```
However, since fc_unzip would exit with error as [***here***](#error) states, I have to follow the compromise solution [***here***](#solution). The following bash commands describe how this compromise solution would be reached:
```bash
# In order to save disk space, extract only .fa files out and delete the intermediate folders
mv 2-asm-falcon/*.fa .
mv 3-unzip/*.fa .
find -maxdepth 1|grep \.\/[0-9].|while read dir;do rm -rf $dir;done
```

<a name="quabus"></a>
## Assembly evaluation by QUAST & BUSCO
```bash
conda activate py2
conda activate busco
```
