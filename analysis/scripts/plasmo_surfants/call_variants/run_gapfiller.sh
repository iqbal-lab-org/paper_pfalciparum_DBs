usage(){
    echo "usage: $0 input_bed input_ref input_bam seed_ins seed_var sample_name"; exit 1
}
input_bed=$1
input_ref=$2
input_bam=$3
seed_ins=$4
seed_var=$5
sample_name=$6

if [[ $# != 6 ]]; then usage; fi

IFS=$'\n';for gene_line in $(cat ${input_bed})
do
    IFS=$'\t'; elems=($gene_line)
    adjusted_start=$((${elems[1]} + 1))
    reg="${elems[0]}:${adjusted_start}-${elems[2]}"
    gene_name=${elems[3]}
    samtools view -h ${input_bam} $reg | samtools sort -n > gene_reads.bam
    GapFiller --query-sam ${input_bam} --seed-sam gene_reads.bam --seed-ins ${seed_ins} --seed-var ${seed_var} --output-prefix ${gene_name}_contigs
    seqkit replace -p "contig" -r "${gene_name}_contig" ${gene_name}_contigs.fasta >> all_gene_contigs.fasta
done
minimap2 -a -c --cs -R "@RG\\tID:1\\tSM:${sample_name}" ${input_ref} all_gene_contigs.fasta | samtools sort > mapped_contigs.bam
samtools index mapped_contigs.bam
bgzip -dc ${input_ref} > input_ref.fa
samtools faidx input_ref.fa
octopus -I mapped_contigs.bam -R input_ref.fa --organism-ploidy 1 --assemble-all -o gapfiller.vcf
