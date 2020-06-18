<h1 align="center">TWO GENOME ASSEMBLY PROJECT</h1>
<p align="left">Litter decomposing fungi represent important agents of soil carbon cycling across terrestrial ecosystems ranging from tundra to tropical forests. The reasons behind the enormous success of litter decomposers are not well known, but it is reasonable to consider that this could be partly due to the diverse plant cell wall decomposition systems these fungi harbor. To examine the diversity of these systems in litter decomposers, we sequenced nine species with diverse habitat preferences from Agaricales, which represents one of the most species-rich order in mushroom forming fungi and in which many litter decomposers are found. All genomes were sequenced using PacBio and the assemblies of seven of these genomes were of good quality. However, the genome assemblies of Mycetinis scorodonius(Mysco, pb_320-2) and Leucopaxillus gentianeus(Leuge, pb_279) did not pass the strict quality standards, that we have set up in this project. The reasons for this could be related to the larger than expected genome size for both species in combination to their dikaryotic state. The successful genome assembly for the two species is of particular interest since both species are found exclusively on conifer litter and therefore, may harbor adaptations related to this recalcitrant type of substrate. Genome assembly tools for dikaryotic genomes are continuously improving, which poses a great opportunity to attempt again the assembly of the two remaining genomes using recently developed tools. </p>

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
***

<a name="intro1"></a>
# Introduction to software which would be used or tried for assembling genomes
Because of the dikaryotic nature of these two fungi, a chronological combination of falcon & falcon-unzip & falcon-phase from  [***pb-assembly***](#pbasm) is supposed to be the most appropriate combo here. Nonetheless, several other popular genome assemblers & polishers ([***miniasm***](#miniasm),[***minipolish***](#minipolish),[***raven***](#raven),[***flye***](#flye),[***canu***](#canu)) which are available for PacBio subreads would be listed and tried as well. (However, [***CANU***](#canu) assembler was ultimately discarded since it had already taken too much time and still hadn't given final assemblies).

Since I have encountered the same issue when using falcon-unzip as the unsolved issue here shows, 

- [Error at 4-polish](https://github.com/PacificBiosciences/FALCON_unzip/issues/159)

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
conda create --name py3 python=3.7
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

<a name="quast"></a>
## QUAST

<a name="bandage"></a>
## Bandage

<a name="busco"></a>
## busco

<a name="pipelines"></a>
# Pipelines of analysis
