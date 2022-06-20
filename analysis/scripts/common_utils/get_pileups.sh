set -eu
BASE_DIR="/hps/nobackup/iqbal/bletcher/plasmo_surfants"

sample=$1
mode=$2
tool=$3
modes=("eval_varcalls" "generational_samples")
choices_ev=("joint_geno" "baseline" "pf6")
choices_gs=("joint_geno" "gapfiller")
usage(){
    echo -e "usage: $0 sample_name mode tool_name;" 
    echo "mode choice: ${modes[@]}"
    echo -e "tool_name choices:"
    echo -e "\t ${choices_ev[@]} (mode = ${modes[0]})"
    echo -e "\t ${choices_gs[@]} (mode = ${modes[1]})"
    exit 0
}
JOINT_GENO_TARGET="gram_jointgeno_ebf7bcd5__pf6__analysis_set_fws95__gapfiller__7__13"
case "$mode" in
    "${modes[0]}")
        case "$tool" in 
            "${choices_ev[0]}")
                TARGET="$JOINT_GENO_TARGET"
                ;;
            "${choices_ev[1]}")
                TARGET=$tool
                ;;
            "${choices_ev[2]}")
                TARGET=$tool
                ;;
            *)
                usage
        esac
        TARGET_DIR="${BASE_DIR}/analysis/outputs/eval_varcalls/induced_refs/pf6/pf6_26_genes/${TARGET}"
        BED_FNAME="${sample}_ir_stats.translated_bed.txt"
        INDUCED_REF="${TARGET_DIR}/${sample}/induced_ref.fa.gz"
        MAPPED_BAM="${TARGET_DIR}/${sample}/mapped.bam"
        ;;
    "${modes[1]}")
        case "$tool" in 
            "${choices_gs[0]}")
                TARGET="$JOINT_GENO_TARGET"
                ;;
            "${choices_gs[1]}")
                TARGET=$tool
                ;;
            *)
                usage
                ;;
        esac
        TARGET_DIR=$(find "${BASE_DIR}/analysis/outputs/generational_samples/${TARGET}" -type d -name "$sample")
        BED_FNAME="pf6_26_genes_induced.bed"
        INDUCED_REF="${TARGET_DIR}/induced_ref.fa.gz"
        MAPPED_BAM="${TARGET_DIR}/induced_ref_mapped.bam"
        ;;
    *)
        usage
        ;;
    esac

BASE_OUTDIR="${BASE_DIR}/tmp_work/induced_ref_mappings"
OUTDIR="${BASE_OUTDIR}/${sample}/${tool}/"
#for sample in PT0177-C PC0302-C PF0304-C QE0438-C PV0287-C PV0245-Cx2 PD1427-C PE0540-C PT0126-C PA0376-C PA0516-C PF0047-C PH0347-C QV0089-C PA0240-C PR0160-C PF0371-C PH1378-C QP0229-C QP0222-C; do
mkdir -p $OUTDIR
induced_ref_out="${OUTDIR}/${sample}_induced_ref.fa.gz"
regions_bed="${TARGET_DIR}/${BED_FNAME}"
regions_pileup_bed="${OUTDIR}/${sample}_pileup_region.bed"
grep -v "_DBL" "${regions_bed}" | awk '{print $1"\t"$2-5000"\t"$3+5000"\t"$4"\t"$5}' > "${regions_pileup_bed}"
cp "${regions_bed}" ${OUTDIR}/${sample}_gene_regions.bed
cp  "${INDUCED_REF}" "$induced_ref_out"
samtools faidx "${induced_ref_out}"
OUTPUT_BAM="${OUTDIR}/${sample}_mapped_to_induced.bam"
samtools view -h -OBAM -L "$regions_pileup_bed" "${MAPPED_BAM}" > ${OUTPUT_BAM}
samtools index ${OUTPUT_BAM}
if [[ -e "${TARGET_DIR}/mapped_contigs.bam" ]]; then
    cp "${TARGET_DIR}"/mapped_contigs.bam* "${OUTDIR}"
fi
#done
