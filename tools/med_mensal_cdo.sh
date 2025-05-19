for file in *.grb2; do
  base_name=$(basename "$file" .grb2)

  cdo timmean "$file" "${base_name}_mean.grb2"
done

cdo mergetime *_mean.grb2 merge_all_files.grb2
