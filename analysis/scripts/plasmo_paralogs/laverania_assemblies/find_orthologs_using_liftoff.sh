set -eu
[[ $# != 4 ]] && echo "usage: $0 gff_file ref_genome_fasta assembly_dir target_gene" \
	&& echo -e "gff_file: features to liftover" \
    && echo -e "assembly_dir: directory containing assembly fastas (.fa.gz)" \
	&& echo -e "target_gene: gene to look for and reconstitute, using 'liftoff', in assemblies in 'assembly_dir'" \
	&& exit 0

gff_file=$1
ref_genome_fasta=$2
assembly_dir=$3
target_gene=$4


mkdir -p tmp
result="${target_gene}_liftoff.fasta"
> ${result}
for asm_name in ${assembly_dir}/*.fa.gz; do
    ID=$(basename $asm_name)
    ID=${ID/.fa.gz/}
    ref="${ID}.fa"
    lifted="tmp/${ID}_liftoff.gff3"
    gzip -dc $asm_name > ${ref}
    liftoff -g ${gff_file} ${ref} ${ref_genome_fasta} -o ${lifted}
    grep "${target_gene};" ${lifted} > lifted.gff3
    bedtools getfasta -fi ${ref} -bed lifted.gff3 | sed 's/>/>'${target_gene}_${ID}' /' >> ${result}
    rm ${ref}*
done
# rm -r tmp
