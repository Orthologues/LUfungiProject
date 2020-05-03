# spades

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
  
