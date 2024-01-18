#!/usr/bin/env bash
# Usage: ./make-methylation_pipeline_manifest.sh 14455_p1

sample=$1
ASD=/net/eichler/vol28/projects/autism_genome_assembly/nobackups

fn=$(echo $sample | cut -f1 -d'_')
sex=$(grep $sample $ASD/sample_info.tab | cut -f2)
mo_sr=/net/eichler/vol28/projects/autism_genome_assembly/nobackups/data/Illumina/WGS/fastq/fofn/${fn}_mo_KG.fastq.fofn
fa_sr=/net/eichler/vol28/projects/autism_genome_assembly/nobackups/data/Illumina/WGS/fastq/fofn/${fn}_fa_KG.fastq.fofn
fq=fastq_fofn/${sample}_fastq.fofn
bm=unmapped_bam_fofn/${sample}_unmapped_bam.fofn

if [ $sex == "M" ]
then
  ref=GRCh38
else
  ref=GRCh38-noY
fi

if echo $sample | grep -qE "_mo|_fa"
then
  echo -e "${sample}\t\t\t\t${ref}\t${fq}\t${bm}\t${fn}"
else
  echo -e "${sample}\t${mo_sr}\t${fa_sr}\t\t${ref}\t${fq}\t${bm}\t${fn}"
fi