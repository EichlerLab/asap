# ASAP (Autism Susceptibility Analysis Pipeline with focus in Structural Variants)
A place to log tasks involved for this project- either executed sequentially or asynchronously.  

##### Table of Contents
* [Inputs](#inputs)
* [QC](#qc)
* [Genome Assembly](#genome-assembly)
* [Alignment](#alignment)
* [Variant calling](#variant-calling)
* [Methylation](#methylation)

## Inputs
### Sample origin/cohort

<div style="text-align: right"> [:arrow_double_up:](#table-of-contents) </div>

## QC
### Fastq.gz input only
* back-reference-qc: use this pipeline for non-human contamination of reads.
### Fasta.gz + its own Illumina input
* merqury: use this tool/pipeline to assess quality of genome assembly
### Bam input
* verifybamid: use this tool/pipeline to assess contamination of non-humanness as well as inter-sample contamination.

<div style="text-align: right"> [:arrow_double_up:](#table-of-contents) </div>
## Genome Assembly
This step produces a fasta file.
* hifiasm: use this pipeline/tool to assemble sample genome
  * **trio-phased requires** parental Illumina data as input
* Version used for all our samples ATM: hifiasm 0.16.1 with just HiFi data.

<div style="text-align: right"> [:arrow_double_up:](#table-of-contents) </div>
## Alignment
This step is produces a BAM.
* HiFi fastq.gz input
  * pbmm2

<div style="text-align: right"> [:arrow_double_up:](#table-of-contents) </div>
## Variant calling
* pbsv: use this with pbmm2 output

<div style="text-align: right"> [:arrow_double_up:](#table-of-contents) </div>
## Methylation
* HiFi data
  * [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools)
* ONT data
  * [modkit](https://github.com/nanoporetech/modkit)

<div style="text-align: right"> [:arrow_double_up:](#table-of-contents) </div>
