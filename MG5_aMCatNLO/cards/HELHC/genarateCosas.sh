# bb_decay="_x0bb"
bb_decay=""
in_folder="ttx0${bb_decay}_build_folder"
folder_strt="mg_pp_ttx0_cosa_"
folder_end="_0123j${bb_decay}_5f"

for i in `seq -1 0.1 1` ; do 
  cp -r ${in_folder} ${folder_strt}${i/.}${folder_end} ; 
  sed -i.Apagar "s/MyCosa/$i/" ${folder_strt}${i/.}${folder_end}/customizecards.dat ; 
done
