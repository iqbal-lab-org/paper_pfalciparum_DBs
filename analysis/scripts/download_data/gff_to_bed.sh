set -eu
[[ $# != 3 ]] && echo "usage: $0 genes_file gff_file output_bed" \
	&& echo -e "genes_file: one gene name per line" \
	&& echo -e "gff_file: genome annotation file" \
	&& exit 0

genes_file=$1
gff_file=$(realpath $2)
output_bed=$3


tmp_file="${output_bed}.tmp"
> $tmp_file
IFS=$'\n'
for line in $(tail -n+2 $genes_file)
	do 
        IFS=$'\t';row=($line)
        gene_name=${row[0]}
        gene_ID=${row[1]}

		found_annot=$(gzip -dc $gff_file | grep -e "ID=$gene_ID;")
		if [[ $(echo $found_annot | wc -l) != 1 ]]; then
			echo "Error: found no, or more than one, annotation line for gene $gene"
			exit 1
		fi
        # gff is [1-based,1-based], and bed is [0-based,0-based); therefore the start pos
        # in gff needs to be decremented by 1.
        echo $found_annot | awk '{$4-=1;print $1"\t"$4"\t"$5"\t""'"$gene_name"'""\t""'"$gene_ID"'"}' >> $tmp_file
done
bedtools sort -i "$tmp_file" > "$output_bed"
rm "$tmp_file"
