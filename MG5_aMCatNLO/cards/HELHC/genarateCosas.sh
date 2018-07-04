in_folder="ttx0_x0bb_build_folder"
folder_strt="mg_pp_ttx0_cosa_"
folder_end="_0123j_5f_x0bb"

for i in `seq -1 0.1 1` ; do 
  cp -r ${in_folder} ${folder_strt}${i/.}${folder_end} ; 
  sed -i.Apagar "s/MyCosa/$i/" ${folder_strt}${i/.}${folder_end}/customizecards.dat ; 
done
