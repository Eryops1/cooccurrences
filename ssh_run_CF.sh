#!/usr/bin/env bash
set -euo pipefail

############################
# Parameters
############################

atlas=17
which_CF="CF_IM_env"

nnull=5000 # for CF_IM_env, CZ+NY: 500. higher for EU cause thats an expensive step
keep_tol=0.05
target_n=1
batch_size=5000
max_batches=10
keep_tol_elev=0.1
cores=4

############################
# Logging
############################

LOGDIR="$HOME/logs/counterfactuals2"
mkdir -p "$LOGDIR"

timestamp=$(date +"%Y%m%d_%H%M%S")
LOGFILE="$LOGDIR/atlas${atlas}_${which_CF}_${timestamp}.log"


############################
# Run
############################

{
    echo "========================================"
    echo "Started : $(date)"
    echo "Host    : $(hostname)"
    echo "Atlas   : $atlas"
    echo "Method  : $which_CF"
    echo "Cores   : $cores"
    echo "========================================"

    export R_LIBS_USER=""

    /home/tietje/mamba/envs/r453/bin/Rscript \
        01_counterfactuals_server.R \
        "$atlas" \
        "$which_CF" \
        "$nnull" \
        "$keep_tol" \
        "$target_n" \
        "$batch_size" \
        "$max_batches" \
        "$cores" \
        "$keep_tol_elev"

    echo
    echo "Finished: $(date)"
} >"$LOGFILE" 2>&1
