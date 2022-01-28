usage (){
echo "usage: $0 bed_fname vcf_fname ref_genome_fname outdir_fname sample_name [1|0] [1|0]"
echo "Takes a bed, a vcf, and a reference genome, and outputs each sequence from each feature in the bed with variants applied, to files in 'outdir_fname'. sample_name gets placed in output fasta headers."
echo "The penultimate param is 1 for placing sample_name in the header of the output sequences or 0 for not. Default 1"
echo "The ultimate param is 1 for applying all calls in the vcf, ignoring sample_name. Useful for template vcf with sample_name=sample. Default 0."
exit 1
}

if [ $# != "5" -a "$#" != "6" -a "$#" != "7" ]; then usage ;fi

bed_fname="$1"
vcf_fname="$2"
ref_genome_fname="$3"
outdir_fname="$4"
sample_name="$5"
sname_in_output="${6:-1}"
apply_all_samples="${7:-0}"

mkdir -p "$outdir_fname"
cp "$vcf_fname" input.vcf.gz
bcftools index input.vcf.gz
bcftools_sample="-s $sample_name"
if [ $apply_all_samples -eq 1 ]; then bcftools_sample="";fi
IFS=$'\n'; for gene_line in $(cat "$bed_fname")
do
    IFS=$'\t'; elems=($gene_line)    
    adj_start=$((${elems[1]} + 1))
    reg="${elems[0]}:${adj_start}-${elems[2]}"
    gene_name=${elems[3]}
    fout=${outdir_fname}/${gene_name}.fa
    new_name="$gene_name"
    if [ $sname_in_output -eq 1 ]; then new_name="${new_name}_${sample_name}"; fi
    samtools faidx "$ref_genome_fname" $reg | bcftools consensus input.vcf.gz $bcftools_sample |
        sed 's/>.*/>'"$new_name"'/' > $fout
done
