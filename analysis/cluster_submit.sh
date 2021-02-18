#!/usr/bin/env bash

set -e

WORKFLOW=$1

function usage(){
    echo "usage: $0 workflow_name"
    exit 1
}

if [[ -z ${WORKFLOW} ]]; then usage; fi

if [[ -z ${TMPDIR} ]]; then TMPDIR="/hps/nobackup/research/iqbal/bletcher/tmp"; fi

# Reduced set of mount points to speed up singularity execution
# mafft requires a writable /scratch directory
SINGULARITY_BINDS="/hps/nobackup/research/iqbal,/nfs/research1/zi,${TMPDIR}:/scratch" # noah
#SINGULARITY_BINDS="/hps/nobackup2/iqbal,/nfs/leia/research/iqbal,${TMPDIR}:/scratch" # yoda
SINGULARITY_ARGS="--contain --bind $SINGULARITY_BINDS"

# In case snakemake --singularity-args fails, can use environment variables too
export SINGULARITY_CONTAIN=TRUE
export SINGULARITY_BINDPATH="$SINGULARITY_BINDS"

LOG_DIR="analysis/logs"
MEMORY=5000

bsub -R "select[mem>$MEMORY] rusage[mem=$MEMORY] span[hosts=1]" \
    -M "$MEMORY" \
    -o "${LOG_DIR}/${WORKFLOW}.o" \
    -e "${LOG_DIR}/${WORKFLOW}.e" \
    -J "${WORKFLOW}_snakemake" \
    snakemake -s analysis/workflows/${WORKFLOW}/Snakefile \
    --profile lsf --verbose --latency-wait 25 \
    --use-singularity #--singularity-args "$SINGULARITY_ARGS"

