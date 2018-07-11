if [[ $# -eq 0 ]]; then
  read -ep "cosa=" cosaVall
else
  cosaVall=$1 
fi
echo "> cosaVall=${cosaVall} "
echo "./submit_gridpack_generation.sh 15000 15000 2nw mg_pp_ttx0_cosa_${cosaVall}_0123j_5f cards/HELHC/mg_pp_ttx0_cosa_${cosaVall}_0123j_5f 2nw\n"
read -p "Press enter to proced!!"
./submit_gridpack_generation.sh 15000 15000 2nw mg_pp_ttx0_cosa_${cosaVall}_0123j_5f cards/HELHC/mg_pp_ttx0_cosa_${cosaVall}_0123j_5f 2nw
