#!/usr/bin/env bash

images_dir=$1
output_dir=$2
parallel=$3
images=`find ${images_dir} -maxdepth 1 -type f`
palette_file=$4

paletteIDs=(`cat ${palette_file} | jq -rc '.[] | .paletteID'`)
rgbLists=(`cat ${palette_file} | jq -rc '.[] | .rgbList'`)

total=${#paletteIDs[*]}

idx=0
while [ ${idx} -lt ${parallel} ]; do
    echo "#!/usr/bin/env bash" > parallel-${idx}.sh
    chmod +x parallel-${idx}.sh
    ((idx++))
done

idx=0
#pairs_num=${#paletteIDs[@]}
IFS=$'\n'
for image in ${images}; do
    file_type=${image##*.}
    group_id=`basename ${image} .${file_type}`
    redirect_no=$((idx % parallel))
    for (( i=0; i<=$(( $total - 1 )); ++i ));
    do
      paletteID=${paletteIDs[i]}
      rgbList=${rgbLists[i]}
      output_path="${output_dir}${group_id}/${group_id}_${paletteID}.${file_type}"
      echo mkdir -p \"`dirname ${output_path}`\" >> parallel-${redirect_no}.sh
      echo ./recolor \"${image}\" \"${output_path}\" \"$rgbList\" 0 >> parallel-${redirect_no}.sh
    done
    ((idx++))
done
