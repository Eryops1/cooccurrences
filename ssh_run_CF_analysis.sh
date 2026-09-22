#!/usr/bin/env bash
set -euo pipefail

############################
# Parameters
############################

atlas_list=(17 26)
which_CF="CF_IM_elev"    # CF_IM, CF_IM_elev, CF_IM_env, clim
ncores=10               # take it easy! starts 4 processes each with N cores, one result file per process
blas_thred=1
omp_thred=1

############################
# Logging
############################

LOGDIR="$HOME/logs/pairwise_processing/"
mkdir -p "$LOGDIR"

############################
# Run
############################

for atlas in "${atlas_list[@]}"; do

    timestamp=$(date +"%Y%m%d_%H%M%S")
    LOGFILE="$LOGDIR/analysis_atlas${atlas}_${which_CF}_${timestamp}.log"

    {
        echo "========================================"
        echo "Started : $(date)"
        echo "Host    : $(hostname)"
        echo "Atlas   : $atlas"
        echo "Method  : $which_CF"
        echo "Cores   : $ncores"
        echo "========================================"

        export R_LIBS_USER=""

        /home/tietje/mamba/envs/r453/bin/Rscript \
            02_counterfactuals_pairwise_analysis_server.R \
            "$atlas" \
            "$which_CF" \
            "$ncores" \
            "$blas_thred" \
            "$omp_thred"

        echo
        echo "Finished: $(date)"
    } >"$LOGFILE" 2>&1 &

done

wait
