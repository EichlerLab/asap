#!/usr/bin/env bash
# Usage: ./get_ont_stats.sh

source ~/.bash_profile
module load miniconda/4.12.0

LRA=/net/eichler/vol28/projects/long_read_archive/nobackups
SCRIPT_PATH=/net/eichler/vol28/software/pipelines/compteam_tools/get_ont_stats.py

for cohort in pop clinical nhp; do
  declare -a sample_array=($(ls -d ${LRA}/${cohort}/*/raw_data/nanopore/*/fastq))

  for str in "${sample_array[@]}"; do
    sample=$(echo "$str" | sed -r "s/(.*${cohort}\/)(.*)(\/raw_data.*)/\2/")
    $SCRIPT_PATH --sample "${sample}" --cohort ${cohort} --prefix ${LRA}
  done

done
