usage (){
echo "usage: $0 bed_fname vcf_fname ref_genome_fname outdir_fname sample_name"
echo "Takes a bed, a vcf, and a reference genome, and outputs each sequence from each feature in the bed with variants applied, to files in 'outdir_fname'. sample_name gets placed in output fasta headers."
exit 1
}

if [[ "$#" -ne 5 ]]; then usage ;fi

bed_fname="$1"
vcf_fname="$2"
ref_genome_fname="$3"
outdir_fname="$4"
sample_name="$5"


mkdir -p "$outdir_fname"
cp "$vcf_fname" input.vcf.gz
bcftools index input.vcf.gz
IFS=$'\n'; for gene_line in $(cat "$bed_fname")
do
    IFS=$'\t'; elems=($gene_line)    
    adj_start=$((${elems[1]} + 1))
    reg="${elems[0]}:${adj_start}-${elems[2]}"
    gene_name=${elems[3]}
    fout=${outdir_fname}/${gene_name}.fa
    samtools faidx "$ref_genome_fname" $reg | bcftools consensus input.vcf.gz |
        sed 's/>.*/>'"${gene_name}_${sample_name}"'/' > $fout
done
