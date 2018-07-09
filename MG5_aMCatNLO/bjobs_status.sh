readarray bjobsTEMP < <(bjobs)

main_job=( $( echo "${bjobsTEMP[@]}" | grep lxplux* | awk '{print $6}' ) )

current_date=`date '+%Y-%m-%d %H:%M:%S'`

echo -e "\n - - - - - - - - - - $current_date - - - - - - - - - - \n"

for i in "${main_job[@]}" ; do 
  echo -e "Main job:\t" $( echo "${bjobsTEMP[@]}" | grep lxplux* | grep $i )
  echo -e "run:\t" $(echo "${bjobsTEMP[@]}" | grep $i | grep "RUN" | wc | awk '{print $1}')
  echo -e "pend:\t" $(echo "${bjobsTEMP[@]}" | grep $i | grep "PEND" | wc | awk '{print $1}')
  echo ""
done

echo -e "\n Total Main:\t" $(echo ${main_job[@]} | wc | awk '{print $2}')
echo -e "\n Total RUN:\t"  $(echo "${bjobsTEMP[@]}" | grep "RUN" | wc | awk '{print $1}')
echo -e "\n Total PEND:\t" $(echo "${bjobsTEMP[@]}" | grep "PEND" | wc | awk '{print $1}')
echo -e "\n Total PEND 1nw:\t" $(echo "${bjobsTEMP[@]}" | grep "PEND" | grep "1nw" | wc | awk '{print $1}')
echo -e "\n Total PEND 2nw:\t" $(echo "${bjobsTEMP[@]}" | grep "PEND" | grep "2nw" | wc | awk '{print $1}')
