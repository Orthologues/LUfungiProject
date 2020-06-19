helpFunction()
{
   echo "--------------------------------Help-Info--------------------------------"
   echo "An integrated pipeline shell script for assembling polished consensus genomes of a diploid species from its pacbio subreads and subsequently do quast analysis for these assemblies"
   echo "Usage: nohup sh ./integrated_falcon_quast.sh -r repository -n index [options] &"
   echo -e "\t-r [STR] Specifies the parent directory of your 'pb-assembly' directory"
   echo -e "\t-n [STR] Specifies the index of your species in PacBio sequencing. For example, if you input '279', its corresponding string would be 'pb_279'"
   echo -e "\t [Options]"
   echo -e "\t-j [INT] Specifies the number of threads used for alignment in pbmm2. 0 means autodetection. (Default = 0)"
   echo -e "\t-k [INT] Specifies the number of threads used for sorting of pbmm2-aligned subreads.0 means 25% of -j, with a maximum of 8.(Default = 0)"     
   echo -e "\t-m [INT] Specifies the memory per thread for sorting. (Default = 768M)"
   echo -e "\t-t [INT] Specifies the number of threads used for variantCaller --algorithm=arrow (Default = 8)"   
   echo "Disclaimer: You must be under a conda environment which has installed pb-assembly(py2 version), pbmm2, pbindex, samtools, genomicconsensus and quast!"
   exit 1
}

while getopts "r:n:j:k:m:t:" opt
do
   case "$opt" in 
      r ) parR="$OPTARG" ;;
      n ) parN="$OPTARG" ;;
      j ) parJ="$OPTARG" ;;
      k ) parK="$OPTARG" ;;
      m ) parM="$OPTARG" ;;
      t ) parT="$OPTARG" ;;
   esac
done

if [[ -z "$parR" ]] || [[ -z "$parN" ]];then
  echo "Some of the compulsory parameters are empty"
  helpFunction
else
  if [[ ! $parR == */ ]];then
     parR+="/"
  fi
  if [[ -z "$parJ" ]];then
     parJ=0
  fi
  if [[ -z "$parK" ]];then
     parK=0
  fi 
  if [[ -z "$parM" ]];then
     parM=768M
  fi
  if [[ -z "$parT" ]];then
     parT=8
  fi
fi
echo "Your input is successful!"
rm ~/${parR}countDone${parN}.txt
touch ~/${parR}countDone${parN}.txt
cd ~/${parR}pb-assembly/
versions=$(find -maxdepth 1 -name "pb_${parN}_v*"|wc -l);
for ((i=1;i<=$versions;i++))
do
( cd ~/${parR}pb-assembly/pb_${parN}_v${i}
  nohup fc_run ../mycfgs/fc_pb_${parN}_v${i}.cfg &> run0.log &
  wait
  nohup fc_unzip.py ../mycfgs/fc_unzip_pb_${parN}.cfg &> run1.std &
  wait
  cat 2-asm-falcon/p_ctg.fa 2-asm-falcon/a_ctg.fa > pb_${parN}_falcon_step2_v${i}.fasta
  cat 3-unzip/all_p_ctg.fa 3-unzip/all_h_ctg.fa >pb_${parN}_falcon_step3_v${i}.fasta
  mv 2-asm-falcon/*.fa .
  mv 3-unzip/*.fa .
  find -maxdepth 1|grep \.\/[0-9].|while read dir;do rm -rf $dir;done
  ls ../../bamfiles/pb_${parN}/*subreads.bam|while read bam;do echo $bam >> pb_${parN}_bam.fofn;done
  nohup pbmm2 align pb_${parN}_falcon_step3_v${i}.fasta pb_${parN}_bam.fofn pb_${parN}_v${i}_aligned.bam --sort -j ${parJ} -J ${parK} -m ${parM} --preset SUBREAD & 
  nohup samtools faidx pb_${parN}_falcon_step3_v${i}.fasta -o pb_${parN}_falcon_step3_v${i}.fasta.fai &
  wait
  nohup pbindex pb_${parN}_v${i}_aligned.bam &
  wait
  nohup arrow pb_${parN}_v${i}_aligned.bam -r pb_${parN}_falcon_step3_v${i}.fasta -o pb_${parN}_step3_v${i}_polished.fastq -j ${parT} --diploid &
  wait
  nohup cat pb_${parN}_step3_v${i}_polished.fastq|paste - - - -|sed 's/^@/>/'|awk '{print $1"\n"$2}' > pb_${parN}_step3_v${i}_polished.fasta &  
  wait
  echo "Done" >> ../../countDone${parN}.txt ) & 
done 
cd ~/${parR}
count=$(cat countDone${parN}.txt|wc -l)
while [ ! "$count" == "$versions" ]
do 
  sleep 10
  count=$(cat countDone${parN}.txt|wc -l)
done
rm ~/${parR}countDone${parN}.txt
touch ~/${parR}countDone${parN}.txt
cd ~/${parR}pb-assembly/
for ((i=1;i<=$versions;i++))
do
( cd ~/${parR}pb-assembly/pb_${parN}_v${i}
  nohup quast.py -o step2_v${i}_quast/  pb_${parN}_falcon_step2_v${i}.fasta -r ../../OriginalAssemblies/pb_${parN}_Leuge.fasta -t 20 &
  nohup quast.py -o step3_v${i}_quast/  pb_${parN}_falcon_step3_v${i}.fasta -r ../../OriginalAssemblies/pb_${parN}_Leuge.fasta -t 20 &
  nohup quast.py -o step4_v${i}_quast/  pb_${parN}_step3_v${i}_polished.fasta -r ../../OriginalAssemblies/pb_${parN}_Leuge.fasta -t 20 & 
  wait
  echo "Done" >> ../../countDone${parN}.txt ) &
done 
cd ~/${parR}
count=$(cat countDone${parN}.txt|wc -l)
while [ ! "$count" == "$versions" ]
do 
  sleep 10
  count=$(cat countDone${parN}.txt|wc -l)
done
cd ~/${parR}pb-assembly/
find -maxdepth 3 -name "report.pdf"|while read pdf;do dir=$(echo $pdf|cut -d / -f 1-3|sed -r 's/quast/quast\//');newname=$(echo $dir|cut -d / -f 2-3|sed -r 's/v[0-9]_//'|tr / _);mv "$pdf" "$dir${newname}.pdf";done
find -maxdepth 3 -name "*quast.pdf" #to check if pdf names are changed properly
