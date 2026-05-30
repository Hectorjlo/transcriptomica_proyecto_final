for file in results/star/DESeq2/*-*.txt; do
name=$(basename $file)
output_path="${file%/$name}"
output_name="${name%.txt}"
output_name="${output_name%%0.[0-9]}_Ensembl_IDs.txt"

awk 'NR > 1 {sub(/\.+[0-9]/, "", $1); print $1}' $file > $output_path/$output_name
done