<h1 align="center">TWO GENOME ASSEMBLY PROJECT</h1>
<p align="left">Litter decomposing fungi represent important agents of soil carbon cycling across terrestrial ecosystems ranging from tundra to tropical forests. The reasons behind the enormous success of litter decomposers are not well known, but it is reasonable to consider that this could be partly due to the diverse plant cell wall decomposition systems these fungi harbor. To examine the diversity of these systems in litter decomposers, we sequenced nine species with diverse habitat preferences from Agaricales, which represents one of the most species-rich order in mushroom forming fungi and in which many litter decomposers are found. All genomes were sequenced using PacBio and the assemblies of seven of these genomes were of good quality. However, the genome assemblies of Mycetinis scorodonius(Mysco, pb_320-2) and Leucopaxillus gentianeus(Leuge, pb_279) did not pass the strict quality standards, that we have set up in this project. The reasons for this could be related to the larger than expected genome size for both species in combination to their dikaryotic state. The successful genome assembly for the two species is of particular interest since both species are found exclusively on conifer litter and therefore, may harbor adaptations related to this recalcitrant type of substrate. Genome assembly tools for dikaryotic genomes are continuously improving, which poses a great opportunity to attempt again the assembly of the two remaining genomes using recently developed tools. </p>

***
# Introduction to software which would be used or tried for assembling genomes
Because of the dikaryotic nature of these two fungi, falcon & falcon-unzip & falcon-phase is supposed to be the most appropriate combo here. Nonetheless, several other popular genome assemblers for PacBio subreads would be listed and tried as well. (However, CANU assembler was ultimately discarded since it really had taken too much time without results)

## blasr
Blasr is a long read aligner which requires a reference genome for input as well. It would be used during the first step 0-phasing step of falcon-unzip of pb-assembly. 

### Installation and help
```bash
conda install -c bioconda blasr
blasr -h
```
### Official Github Page
- [OfficialLink](https://github.com/PacificBiosciences/blasr)

## pbmm2
pbmm2 is a SMRT C++ wrapper for minimap2's C API. Its purpose is to support native PacBio in- and output, provide sets of recommended parameters, generate sorted output on-the-fly, and postprocess alignments. Sorted output can be used directly for polishing using GenomicConsensus, if BAM has been used as input to pbmm2. Benchmarks show that pbmm2 outperforms BLASR in sequence identity, number of mapped bases, and especially runtime. pbmm2 is the official replacement for BLASR.

### Installation and help
```bash
conda install -c bioconda pbmm2
pbmm2 -h
```
### Official Github Page
- [OfficialLink](https://github.com/PacificBiosciences/pbmm2)

## GenomicConsensus
The GenomicConsensus package provides the variant Caller tool, which allows you to apply the Quiver or Arrow algorithm to mapped PacBio reads to get consensus and variant calls. It will be used as one of the dependencies of falcon-unzip of pb-assembly. It would require unaligned .bam files as raw sequence data with 'signal pulse' information and be used for phased-polishing with Arrow.

### Installation and help

```bash
conda install -c bioconda genomicconsensus python=2
quiver -h
arrow -h
```
### Official Github Page
- [OfficialLink](https://github.com/PacificBiosciences/GenomicConsensus)

## pb-assembly
pb-assembly is an official PacBio combo of three genome assembly tools - falcon, falcon-unzip and falcon-phase. It provides an integrated workflow to assemble diploid genomes from subread files whose formats are .fasta and .bam.
### Installation and help

```bash
conda create -n denovo
conda activate denovo
conda install -c bioconda pb-assembly
fc_run.py -h
```
### Official Github Page
- [OfficialLink](https://github.com/PacificBiosciences/pb-assembly)

## miniasm
Miniasm is a rather particular long-read assembler as it doesn't include a consensus step. The resulting contigs are just merged erroneous long reads and still contains many sequencing errors. Produced contigs are structually correct, but at the nucleotide level, there're many mismatches and indels. Its output format is .gfa.

### Official Github Page
- [OfficialLink](https://github.com/lh3/miniasm)

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
- [OfficialLink](https://github.com/rrwick/Minipolish#installation)

## raven

## flye

## canu

# Introduction to software which would be used or tried for evaluation of raw data & assembly evaluation

## fastqc

## multiqc

## QUAST

## Bandage

## busco

