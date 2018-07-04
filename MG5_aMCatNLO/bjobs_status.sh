readarray bjobsTEMP < <(bjobs)

main_job=$( echo "${bjobsTEMP[@]}" | grep lxplux* | awk '{print $6}' )

for i in $main_job ; do 
  echo -e "Main job:\t" $( echo "${bjobsTEMP[@]}" | grep lxplux* | grep $i )
  echo -e "run:\t" $(echo "${bjobsTEMP[@]}" | grep $i | grep "RUN" | wc | awk '{print $1}')
  echo -e "pend:\t" $(echo "${bjobsTEMP[@]}" | grep $i | grep "PEND" | wc | awk '{print $1}')
  echo ""
done

echo -e "\n Total Main:\t" $(echo ${main_job[@]} | wc | awk '{print $2}')
echo -e "\n Total RUN:\t"  $(echo "${bjobsTEMP[@]}" | grep "RUN" | wc | awk '{print $1}')
echo -e "\n Total PEND:\t" $(echo "${bjobsTEMP[@]}" | grep "PEND" | wc | awk '{print $1}')
