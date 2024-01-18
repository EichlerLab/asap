#!/usr/bin/env bash
# Usage: ./run_hifiasm.sh 30

JOB_COUNT=$1

shift

DRMAA_ARGS=" -V -cwd -j y -o ./log -e ./log -l h_rt={resources.hrs}:00:00 -l mfree={resources.mem}G -pe serial {threads} -w n -S /bin/bash"

# Load the module where snakemake is installed
module load miniconda/4.12.0

# Make log directory
mkdir -p log

# Run
snakemake \
    -s /net/eichler/vol28/software/pipelines/hifiasm-smk/hifi/Snakefile \
    --restart-times 3 \
    --use-envmodules \
    --drmaa "${DRMAA_ARGS}" \
    --jobname "{rulename}.{jobid}" \
    --printshellcmds \
    --keep-going \
    -j "${JOB_COUNT}" \
    "$@"
