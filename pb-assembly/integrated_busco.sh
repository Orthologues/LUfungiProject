#This script automates busco analysis of existent falcon assemblies 
helpFunction()
{
   echo "--------------------------------Help-Info--------------------------------"
   echo "Usage: nohup sh ./integrated_busco.sh -r repository -n index -k Kingdom &"
   echo "An integrated pipeline shell script to do busco analysis for falcon-assemblies of different steps after running integrated_falcon_quast.sh"
   echo -e "\t-r [STR] Specifies the parent directory of your 'pb-assembly' directory"
   echo -e "\t-n [STR] Specifies the index of your species in PacBio sequencing. For example, if you input '279', its corresponding string would be 'pb_279'"
   echo -e "\t-k [STR] Specifies the kingdom name of your busco dataset. For example, 'fungi_odb10'\n"
   echo "Disclaimer: You must have installed conda and have run these commands before using this script:"
   echo -e "\tconda create -n busco -c bioconda -c conda-forge busco=4.0.6 python=3.7"
   echo -e "\tconda deactivate"
   echo -e "\tconda activate busco"
   exit 1
}

while getopts "r:n:k:" opt
do
   case "$opt" in 
      r ) parR="$OPTARG" ;;
      n ) parN="$OPTARG" ;;
      k ) parK="$OPTARG" ;;
   esac
done

if [[ -z "$parR" ]] || [[ -z "$parN" ]] || [[ -z "$parK" ]];then
  echo "Some of the compulsory parameters are empty"
  helpFunction
else
  if [[ ! $parR == */ ]];then
     parR+="/"
  fi
fi

if [[ ! -f ~/${parR}countDone${parN}.txt ]];then
  echo -e "You haven't finished running 'integrated_falcon_quast.sh' yet!"
  helpFunction
else
  cd ~/${parR}pb-assembly/
  count=$(cat ../countDone${parN}.txt|wc -l)
  versions=$(find -maxdepth 1 -name "pb_${parN}_v*"|wc -l);
  if [[ ! "$count" == "$versions" ]];then 
    echo -e "You haven't finished running 'integrated_falcon_quast.sh' yet!"
    helpFunction
  else 
    echo -e "This script has begun to run...\n"
    rm countDone${parN}_busco.txt
    touch countDone${parN}_busco.txt
    for ((i=1;i<=$versions;i++))
    do
    ( cd ~/${parR}pb-assembly/pb_${parN}_v${i}
      nohup busco -m genome -i pb_${parN}_falcon_step2_v${i}.fasta -o step2_busco -l ${parK} & sleep 10
      nohup busco -m genome -i pb_${parN}_falcon_step3_v${i}.fasta -o step3_busco -l ${parK} & sleep 10
      nohup busco -m genome -i pb_${parN}_step3_v${i}_polished.fasta -o step4_busco -l ${parK} & sleep 10
      wait
      echo "Done" >> ../countDone${parN}_busco.txt ) &
    done 
    cd ~/${parR}pb-assembly
    countBusco=$(cat countDone${parN}_busco.txt|wc -l)
    while [ ! "$countBusco" == "$versions" ]
    do 
      sleep 10
      countBusco=$(cat countDone${parN}_busco.txt|wc -l)
    done
    mkdir busco_plots${parN}
    cd busco_plots${parN}
    ls ../*${parN}*/step*_busco/*.txt|while read txt;do name=$(echo $txt|cut -d / -f 2-3|tr -d \/|sed 's/step/_step/');echo $name;mkdir $name;cp $txt $name;done
    nohup sh -c 'ls *busco/|grep ^pb|sed 's/://'|while read dir;do nohup generate_plot.py -wd $dir;done' & wait
  fi  
fi
