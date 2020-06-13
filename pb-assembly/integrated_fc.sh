helpFunction()
{
   echo "-----------------------------Help-Info--------------------------------"
   echo "An integrated pipeline shell script for assembling polished consensus genomes of a diploid species from its pacbio subreads"
   echo "Usage: sh ./integrated_fc.sh -r repository -n index [options] -h &"
   echo -e "\t-r [STR] Specifies the parent directory of your 'pb-assembly' directory"
   echo -e "\t-n [STR] Specifies the index of your species in PacBio sequencing. For example, if you input '279', its corresponding string would be 'pb_279'"
   echo -e "\t [Options]"
   echo -e "\t-j [INT] Specifies the number of threads used for alignment in pbmm2. 0 means autodetection. (Default = 0)"
   echo -e "\t-k [INT] Specifies the number of threads used for sorting of pbmm2-aligned subreads.0 means 25% of -j, with a maximum of 8.(Default = 0)"     
   echo -e "\t-m [INT] Specifies the memory per thread for sorting. (Default = 768M)"
   echo -e "\t-t [INT] Specifies the number of threads used for variantCaller --algorithm=arrow (Default = 8)"     
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