#!/usr/bin/env bash
# Usage: ./run.sh 30

JOB_COUNT=$1

shift

DRMAA_ARGS=" -V -cwd -j y -o ./log -e ./log -l h_rt={resources.hrs}:00:00 -l mfree={resources.mem}G -pe serial {threads} -w n -S /bin/bash"

# Load the module where snakemake is installed
module load miniconda/4.12.0

# Make log directory
mkdir -p log

# Run
snakemake \
    -s qc-sex.smk \
    --configfile config.yaml \
    --restart-times 3 \
    --use-envmodules \
    --jobname "{rulename}.{jobid}" \
    --drmaa "${DRMAA_ARGS}" \
    --printshellcmds \
    --keep-going \
    -w 5 \
    -j "${JOB_COUNT}" \
    "$@"
