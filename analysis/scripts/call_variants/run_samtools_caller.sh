#!/usr/bin/env bash
set -e

args=("$@")
if [ ${#args[@]} -ne 4 ]; then echo "usage: $0 bwa_index_prefix reads_files output_vcf num_threads"; exit 0; fi
idx_prefix="${args[0]}"
reads_files="${args[1]}"
output_vcf="${args[2]}"
threads="${args[3]}"
bwa mem -t $threads $idx_prefix $reads_files > mapped.sam
samtools sort -o mapped.bam -O BAM mapped.sam
samtools index mapped.bam
bcftools mpileup -f $idx_prefix mapped.bam | bcftools call -O z -o $output_vcf -vm --ploidy 1
