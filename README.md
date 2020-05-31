# GenomicConsensus
The GenomicConsensus package provides the variant Caller tool, which allows you to apply the Quiver or Arrow algorithm to mapped PacBio reads to
get consensus and variant calls. It will be used as one of the dependencies of pb-assembly.

## Installation and help

```bash
conda install -c bioconda genomicconsensus python=2
quiver -h
arrow -h
```
## Official Github Page
[OfficialLink](https://github.com/PacificBiosciences/GenomicConsensus)

# miniasm
Miniasm is a rather particular long-read assembler as it doesn't include a consensus step. The resulting contigs are just merged erroneous long reads and still contains many sequencing errors. Produced contigs are structually correct, but at the nucleotide level, there're many mismatches and indels. Its output format is .gfa.

## Official Github Page
[OfficialLink](https://github.com/lh3/miniasm)

# minipolish
Minipolish is a combined tool for .gfa genome assembly polishing, especially for those which were produced by miniasm. It combines racon and minimap and can run for multiple rounds. It supports both ONT and PACBIO subreads as input.

## Installation and help

```bash
conda create --name py3 python=3.7
conda activate py3
conda install -c bioconda minipolish
minipolish -h
```

## Official github page
[OfficialLink](https://github.com/rrwick/Minipolish#installation)

# blasr

# raven

# flye

# canu

# QUAST

# Bandage

# pb-assembly
