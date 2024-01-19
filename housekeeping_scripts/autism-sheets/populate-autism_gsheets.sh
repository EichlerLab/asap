#!/usr/bin/env bash
# Usage: ./populate-autism_gsheets.sh

# Load the required executables
source ~/.bash_profile
module load miniconda/4.12.0

for t in nanopore PacBio_HiFi other
do
  if [ $t == "PacBio_HiFi" ]
  then
    nt="hifi"
  elif [ $t == "nanopore" ]
  then
    nt="ont"
  else
    nt=$t
  fi

  target_file=db-${nt}.tsv.gz
  target_sheet=db-${nt}
  ~/nobackups/scripts/autism-sheets.py prepare\
    --sample_file /net/eichler/vol28/projects/autism_genome_assembly/nobackups/sample_info.tab \
    --tech $t \
    --illumina_fofn_prefix /net/eichler/vol28/projects/autism_genome_assembly/nobackups/data/Illumina/WGS/fastq/fofn \
    --merqury_prefix /net/eichler/vol28/projects/autism_genome_assembly/nobackups/qc/merqury/merqury/ \
    --vbi_prefix /net/eichler/vol28/projects/autism_genome_assembly/nobackups/qc/verifybamid \
    --methyl_prefix /net/eichler/vol28/projects/autism_genome_assembly/nobackups/methylation/new-run/ \
    --brqc_prefix /net/eichler/vol28/projects/autism_genome_assembly/nobackups/qc/back-reference-qc/ \
    --outpath $target_file

  until [ -f $target_file ]
  do
     sleep 1
  done

 ./autism-sheets.py sheets $target_file --sheet_name $target_sheet

done
