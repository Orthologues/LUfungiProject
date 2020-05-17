# SPAdes
Spades is a popular genome assembler. One of its advantages is that it alters 
the k-mer size depending on the coverage of the genome. For low coverage re-
gions  k-mer size should be low and for high coverage regions high. From supplementary material it can be learned that '--pacbio' option should be used

## Installation and testing

```bash
cd ~/bin
wget http://cab.spbu.ru/files/release3.14.1/SPAdes-3.14.1-Linux.tar.gz
tar -xvzf SPAdes-3.14.1-Linux.tar.gz 
cd ~/bin/SPAdes-3.14.1-Linux/bin
spades.py -h
```

## Instruction of options
[OfficialLink](http://cab.spbu.ru/files/release3.14.1/manual.html)
  
# GenomicConsensus

The GenomicConsensus package provides the variant Caller tool, which allows you to apply the Quiver or Arrow algorithm to mapped PacBio reads to
get consensus and variant calls

## Installation and help

```bash
conda install -c bioconda genomicconsensus python=2
quiver -h
arrow -h
```

## Official Github Page
[OfficialLink](https://github.com/PacificBiosciences/GenomicConsensus)

# wtdbg2 
It's a de novo sequence assembler for long noisy reads produced by PacBio or Oxford Nanopore Technologies (ONT) which can be used for genome polishing as well

## Official Github Page
[OfficialLink](https://github.com/ruanjue/wtdbg2)

# racon
It's a Consensus module for raw de novo DNA assembly of long uncorrected reads. Racon can be used as a polishing tool after the assembly with either Illumina data or data produced by third generation of sequencing.

## Official Github Page
[OfficialLink](https://github.com/isovic/racon)

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


