#!/usr/bin/env bash
# Usage: ./populate-gsheets_sequencing-summary.sh

# Load the required executables
source ~/.bash_profile
module load miniconda/4.12.0

for s in {clinical,pop,nhp}:{nanopore,PacBio_HiFi}
do
  c=$(echo $s | cut -f1 -d':')
  t=$(echo $s | cut -f2 -d':')
  if [ $t == "PacBio_HiFi" ]
  then
    nt="hifi"
  else
    nt="ont"
  fi

  target_file=${nt}-${c}.tsv.gz
  target_sheet=${nt}-${c}
  ./prepare.py \
    --tech $t \
    --cohort $c \
    --proj_dir /net/eichler/vol28/projects/long_read_archive/nobackups/ \
    --outpath $target_file


  until [ -f $target_file ]
  do
     sleep 1
  done

 ./sheets.py $target_file --sheet_name $target_sheet

done
