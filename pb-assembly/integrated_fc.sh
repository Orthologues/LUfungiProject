helpFunc()
{
   echo "-------------------------------------------------------------------"
   echo "An integrated pipeline shell script for assembling polished consensus genomes of a diploid species from its pacbio subreads"
   echo "Usage: sh ./integrated_fc.sh -n index [options] &"
   echo -e "\t-n [STR] Specifies the index of your species in PacBio sequencing. For example, if you input '279', its corresponding string would be 'pb_279'"
   echo -e "\t-j [INT] Specifies the number of threads used for alignment in pbmm2. 0 means autodetection. (Default = 0)"
   echo -e "\t-J [INT] Specifies the number of threads used for sorting of pbmm2-aligned subreads.0 means 25% of -j, with a maximum of 8.(Default = 0)"     
   
}
