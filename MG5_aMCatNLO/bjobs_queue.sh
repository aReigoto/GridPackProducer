readarray bjobsTEMP < <(bjobs)

main_job=( $( echo "${bjobsTEMP[@]}" | grep lxplux* | awk '{print $6}' ) )

# for i in "${main_job[@]}" ; do
for i in "${bjobsTEMP[@]}" ; do
  echo "$i" | grep "PEND" | grep "1nw" | awk '{print $1}' | xargs -I {} bswitch 2nw {}
  # echo "${bjobsTEMP[@]}" | grep $i
  # echo "${bjobsTEMP[@]}" | grep $i | grep "PEND" | grep "1nw" |  awk '{print $1}' | xargs -I {} bswitch 2nw {}
  # echo "${bjobsTEMP[@]}" | grep $i | grep "PEND" | grep "1nw" |  awk '{print $1}'
done


