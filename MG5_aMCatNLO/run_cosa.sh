cosaVall=$1 && echo "> cosaVall=${cosaVall} "
read tempVar
./submit_gridpack_generation.sh 15000 15000 2nw mg_pp_ttx0_cosa_${cosaVall}_0123j_5f cards/HELHC/mg_pp_ttx0_cosa_${cosaVall}_0123j_5f 1nw
