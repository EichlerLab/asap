#!/usr/bin/env bash
# Usage: ./sfari-make-metadata.sh ssc_name:sex:family:member path/to/lra

ssc_name=$(echo $1 | cut -f1 -d':')

# Check if directory exists
if [ ! -d "./raw_data/${ssc_name}/" ]
then
  echo "${ssc_name} does not exist- skipping." 1>&2
  exit 1
fi

sex=$(echo $1 | cut -f2 -d':')
family=$(echo $1 | cut -f3 -d':')
member=$(echo $1 | cut -f4 -d':')

lra=$2

# ont
nanopore_runids=$(ls -d ./raw_data/${ssc_name}/ont/STD/bam/* | rev | cut -f1 -d'/' | rev)

# PacBio_HiFi
hifi_runids=$(ls ./raw_data/${ssc_name}/PacBio_HiFi/*.bam | xargs -i basename {} | sed 's/\.bam//g')

#echo -e "technology\trun_identifier\tcell_type\tssc_identifier\tsex\tfamily_id\tfamily_role"
for nano in $nanopore_runids
do
  source="cell-line"
  echo -e "nanopore\t$nano\t${source}\t$ssc_name\t$sex\t$family\t$member"
done

for hifi in $hifi_runids
do
  if ! echo $hifi | grep -q "fail" && ! echo $hifi | grep -qE "[0-9]{1,5}_[0-9]{1,5}.reads"
  then

    source="blood"

    # Identify true source
    module load samtools/1.14

    fp=${lra}/${family}_${member}/raw_data/PacBio_HiFi/${hifi}.bam
    if [ -f $fp ]
    then
      out=$(samtools view -H $fp | grep -Eo "SM\:.*PM")
      source_name=$(echo $out | sed -r 's/SM\://; s/PM//' | tr -d ' ')
      if $(echo $source_name | grep -q "SSC")
      then
        source="cell-line"
      fi
    else
      echo "There is no corresponding bam: ${fp}" 1>&2
      exit 1
    fi

    echo -e "PacBio_HiFi\t$hifi\t$source\t$ssc_name\t$sex\t$family\t$member"
  fi
done
