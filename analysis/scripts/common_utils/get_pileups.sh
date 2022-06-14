sample=$1
tool=$2
choices=("joint_geno" "baseline" "pf6")
usage(){
    echo "usage: $0 sample_name tool_name; \n tool_name choice: ${choices[@]}"
    exit 0
}
case "$tool" in 
    "${choices[0]}")
        TARGET="gram_jointgeno_ebf7bcd5__pf6__analysis_set_fws95__gapfiller__7__13"
        ;;
    "${choices[1]}")
        TARGET=$tool
        ;;
    "${choices[2]}")
        TARGET=$tool
        ;;
    *)
        usage
        ;;
esac

BASE_DIR="/hps/nobackup/iqbal/bletcher/plasmo_surfants"
TARGET_DIR="${BASE_DIR}/analysis/outputs/eval_varcalls/induced_refs/pf6/pf6_26_genes/${TARGET}"
BASE_OUTDIR="${BASE_DIR}/tmp_work/induced_ref_mappings"
OUTDIR="${BASE_OUTDIR}/${sample}/${tool}/"
#for sample in PT0177-C PC0302-C PF0304-C QE0438-C PV0287-C PV0245-Cx2 PD1427-C PE0540-C PT0126-C PA0376-C PA0516-C PF0047-C PH0347-C QV0089-C PA0240-C PR0160-C PF0371-C PH1378-C QP0229-C QP0222-C; do
mkdir -p $OUTDIR
induced_ref="${OUTDIR}/${sample}_induced_ref.fa.gz"
mapped_bam="${TARGET_DIR}/${sample}/mapped.bam"
regions_bed="${TARGET_DIR}/${sample}_ir_stats.translated_bed.txt"
regions_pileup_bed="${OUTDIR}/${sample}_pileup_region.bed"
grep -v "_DBL" "${regions_bed}" | awk '{print $1"\t"$2-5000"\t"$3+5000"\t"$4"\t"$5}' > "${regions_pileup_bed}"
cp "${regions_bed}" ${OUTDIR}/${sample}_gene_regions.bed
cp "${TARGET_DIR}/${sample}/induced_ref.fa.gz" "$induced_ref"
samtools faidx "${induced_ref}"
OUTPUT_BAM="${OUTDIR}/${sample}_mapped_to_induced.bam"
samtools view -h -OBAM -L "$regions_pileup_bed" "${mapped_bam}" > ${OUTPUT_BAM}
samtools index ${OUTPUT_BAM}
#done
