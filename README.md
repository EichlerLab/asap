# ASAP
Autism Susceptibility Analysis Pipeline with focus in Structural Variants.

This repo logs tasks involved for this project- either executed sequentially or asynchronously.  

##### Table of Contents
* [Inputs](#inputs)
* [QC](#qc)
* [Genome Assembly](#genome-assembly)
* [Alignment](#alignment)
* [Variant calling](#variant-calling)
* [Getting a refined call set](#refined-callset)
* [Methylation](#methylation)
* [FAQ](#faq)

## Inputs
### Sample origin/cohort

[:arrow_double_up:](#table-of-contents)
## QC
### Fastq.gz input only
* back-reference-qc: use this pipeline for non-human contamination of reads.
### Fasta.gz + its own Illumina input
* merqury: use this tool/pipeline to assess quality of genome assembly
### Bam input
* verifybamid: use this tool/pipeline to assess contamination of non-humanness as well as inter-sample contamination.

[:arrow_double_up:](#table-of-contents)
## Genome Assembly
This step produces a fasta file.
* hifiasm: use this pipeline/tool to assemble sample genome
  * **trio-phased requires** parental Illumina data as input
* Version used for all our samples ATM (Jan 16, 2024): hifiasm 0.16.1 with just HiFi data.

[:arrow_double_up:](#table-of-contents)
## Alignment
This step is produces a BAM.
* HiFi fastq.gz input
  * pbmm2

[:arrow_double_up:](#table-of-contents)

## Variant calling
* pbsv: use this with pbmm2 output

[:arrow_double_up:](#table-of-contents)

## Refined callset
### SVPOP sampleset merging of SVs
### Discovery curve
### DNM validation

[:arrow_double_up:](#table-of-contents)

## Methylation
This step produces a methylation bed file and
* HiFi data
  * [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools)
* ONT data
  * [modkit](https://github.com/nanoporetech/modkit)

[:arrow_double_up:](#table-of-contents)

## FAQ
1. I think the sample is contaminated, how should I investigate?
   1. Try verifybamid either for each cell *or* for the entire sample aligned against the same reference as its Illumina bam.
   2. Try back-reference-qc to see how much non-human contamination there is and remove it from the fastq.gz
2. I think the sample is of different sex, what to do?
   1. Try looking at coverage over the SRY gene for quick visualization
   2. Estimate the sex using [this script](https://github.com/projectoriented/bio-utils/blob/main/sex-estimator.py) either for each cell or the entire sample. The input must be a BAM.

[:arrow_double_up:](#table-of-contents)