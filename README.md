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
* [Analyses](#analyses)
* [FAQ](#faq)
* [Paths](#paths)
* [Housekeeping](#housekeeping)

## Inputs
### Sample origin/cohort
This pertains to SSC cohort.

| Count | Sex (proband-sibling) | Family type |
|:-----:|:---------------------:|------------:|
|  11   |          F-F          |        quad |
|  17   |          F-M          |        quad |
|   4   |          M-F          |        quad |
|   5   |          M-M          |        quad |
|  10   |           F           |        trio |
|   2   |           M           |        trio |
|   1   |         F-F-F         |       penta |

Separately, there are three SAGE trios: BK143, BK196, BK486.
Sample sheet is here: `/net/eichler/vol28/projects/autism_genome_assembly/nobackups/sample_info.tab`

[:arrow_double_up:](#table-of-contents)
## QC
#### Fastq.gz input only
* back-reference-qc: use this pipeline for non-human contamination of reads.
#### Fasta.gz + its own Illumina input
* merqury: use this tool/pipeline to assess quality of genome assembly
#### Bam input
* verifybamid: use this tool/pipeline to assess contamination of non-humanness as well as inter-sample contamination.

[:arrow_double_up:](#table-of-contents)
## Genome Assembly
This step produces a fasta file.
* hifiasm: use this pipeline/tool to assemble sample genome
  * **trio-phased requires** parental Illumina data as input
* Version used for all our samples ATM (Jan 16, 2024): hifiasm 0.16.1 with just HiFi data.

[:arrow_double_up:](#table-of-contents)
## Alignment
This step is produces a BAM. And can be achieved via: https://eichlerlab.gs.washington.edu/help/wiki/doku.php?id=users:wharvey:align_all
* HiFi fastq.gz input
  * [pbmm2](https://github.com/PacificBiosciences/pbmm2)
* Illumina fastq.gz input

[:arrow_double_up:](#table-of-contents)

## Variant calling
* [pbsv](https://github.com/PacificBiosciences/pbsv): use pbmm2 output as input for this
* [PAV](https://github.com/EichlerLab/pav): use hifiasm assembly output for this
  * [instructions](notes/pav.md)

[:arrow_double_up:](#table-of-contents)

## Refined callset
The steps here are by sequential order.

#### SUBSEQ for SV validation
#### SVPOP sampleset merging of SVs
* [SVPOP](https://github.com/EichlerLab/svpop): [instructions here](notes/svpop.md)
#### Discovery curve
#### DNM validation

[:arrow_double_up:](#table-of-contents)

## Methylation
This step produces a methylation bed file and bigwig files of the beds. The pipeline for both technologies are here: https://eichlerlab.gs.washington.edu/help/wiki/doku.php?id=users:lettucerap:methylation
* HiFi data
  * [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools)
* ONT data
  * [modkit](https://github.com/nanoporetech/modkit)

[:arrow_double_up:](#table-of-contents)

## Analyses
#### Contiguous chromosome X

#### DNM visualization
https://eichlerlab.gs.washington.edu/help/wiki/doku.php?id=users:lettucerap:ava_svbyeye
#### DMR identification

## FAQ
1. I think the sample is contaminated, how should I investigate?
   1. Try verifybamid either for each cell *or* for the entire sample aligned against the same reference as its Illumina bam.
   2. Try back-reference-qc to see how much non-human contamination there is and remove it from the fastq.gz
2. I think the sample is of different sex, what to do?
   1. Try looking at coverage over the SRY gene for quick visualization
   2. Estimate the sex using [this script](https://github.com/projectoriented/bio-utils/blob/main/sex-estimator.py) either for each cell or the entire sample. The input must be a BAM.

[:arrow_double_up:](#table-of-contents)

## Paths
```shell
LRA=/net/eichler/vol28/projects/long_read_archive/nobackups
sample=14455_p1

# File of File Names (fofn)
# hifi
ls -lrtha $LRA/clinical/${sample}/raw_data/PacBio_HiFi/fofn/ccs/fastq.fofn
# ont
PATH="/net/eichler/vol28/software/pipelines/compteam_tools:$PATH"
make_ont_fofn.py --sample ${sample} \
  --proj_dir $LRA/clinical \
  --output ${sample}.fofn \
  --filter_string 'lib=STD;model=sup.*;bc=guppy;ver=6;ftype=bam'

# hifiasm assemblies
ls -lrtha $LRA/clinical/${sample}/assemblies/hifiasm
```

https://eichlerlab.gs.washington.edu/help/wiki/doku.php?id=users:lettucerap:make_ont_fofn

[:arrow_double_up:](#table-of-contents)

## Housekeeping
1. SFARI data deposition
2. Scripts that routinely updates the google sheets.

[:arrow_double_up:](#table-of-contents)
