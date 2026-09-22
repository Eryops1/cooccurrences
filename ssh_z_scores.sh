#!/usr/bin/env bash
set -euo pipefail

############################
# Parameters
############################

ncores=150 # no more than this, uses all CPU
blas_thred=1
omp_thred=1
combs="all"
combs_n=200

############################
# Logging
############################

LOGDIR="$HOME/logs/z_scores"
mkdir -p "$LOGDIR"

timestamp=$(date +"%Y%m%d_%H%M%S")
LOGFILE="$LOGDIR/z_scores_${timestamp}.log"


############################
# Run
############################

{
    echo "========================================"
    echo "Started : $(date)"
    echo "Host    : $(hostname)"
    echo "Cores   : $ncores"
    echo "Pairs   : $combs"
    echo "Pairs N : $combs_n"
    echo "========================================"

    export R_LIBS_USER=""

    /home/tietje/mamba/envs/r453/bin/Rscript \
        03_z_scores_server.R \
        "$ncores" \
        "$blas_thred" \
        "$omp_thred" \
        "$combs" \
        "$combs_n"

    echo
    echo "Finished: $(date)"
} >"$LOGFILE" 2>&1
