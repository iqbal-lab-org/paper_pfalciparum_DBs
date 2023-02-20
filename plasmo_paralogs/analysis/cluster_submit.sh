#!/usr/bin/env bash

set -e

WORKFLOW=$1

function usage(){
    echo "usage: $0 workflow_name [no_singu]"
    exit 1
}

if [[ -z ${WORKFLOW} ]]; then usage; fi

if [[ -z ${TMPDIR} ]]; then TMPDIR="/tmp"; fi

# Reduced set of mount points to speed up singularity execution
# mafft requires a writable /scratch directory
SINGULARITY_BINDS="/hps/nobackup/iqbal,/nfs/research/zi,${TMPDIR}:/scratch,${TMPDIR}:/tmp" # codon
#SINGULARITY_BINDS="/hps/nobackup/research/iqbal,/nfs/research1/zi,${TMPDIR}:/scratch" # noah
#SINGULARITY_BINDS="/hps/nobackup2/iqbal,/nfs/leia/research/iqbal,${TMPDIR}:/scratch" # yoda
SINGULARITY_ARGS="--contain --bind $SINGULARITY_BINDS"

export SINGULARITY_CONTAIN=TRUE
export SINGULARITY_BINDPATH="$SINGULARITY_BINDS"

LOG_DIR="analysis/logs"
mkdir -p "${LOG_DIR}"
MEMORY=5000

cmd='bsub -R "select[mem>$MEMORY] rusage[mem=$MEMORY] span[hosts=1]" \
    -M "$MEMORY" \
    -o "${LOG_DIR}/${WORKFLOW}.o" \
    -e "${LOG_DIR}/${WORKFLOW}.e" \
    -J "${WORKFLOW}_snakemake" \
    snakemake -s analysis/workflows/${WORKFLOW}/Snakefile \
    --profile lsf --verbose --latency-wait 120 --keep-going'
if [[ "$2" != "no_singu" ]]; then cmd="${cmd} --use-singularity"; fi 
#Note: pass args to singu via '--singularity-args "$SINGULARITY_ARGS"'.
# The environment variables (SINGULARITY_CONTAIN, SINGULARITY_BINDPATH) replace need for this.

echo "$cmd"
eval "$cmd"
