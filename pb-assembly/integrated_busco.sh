#This script automates busco analysis of existent falcon assemblies 
helpFunction()
{
   echo "--------------------------------Help-Info--------------------------------"
   echo "Usage: sh ./integrated_busco.sh -r repository -n index -k Kingdom &"
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
helpFunction
