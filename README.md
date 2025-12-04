# ASAP
Autism Susceptibility Analysis Pipeline with a focus on Structural Variants (SVs). This repository documents the tasks involved in this project, which may be executed either sequentially or asynchronously. The approach used for rare variant or pathogenic candidate discovery in this study can be applied broadly to families affected by any rare disease.

##### System Requirements 
Hardware requirements: Any Processor capable of running x86_64 architecture and at least 128GB of memory. Some steps can process samples in parallel, while the steps that handle all samples together scale logarithmically with sample size.
Software requirements: The developed code mainly depends on the Python3 scientific stack and has been tested on the following system: Ubuntu 22.04.

##### Table of Contents
* [Sample](#inputs)
* [QC](#qc)
  * [back-reference-qc](#back-reference-qc)
  * [ntsm](#ntsm)
  * [VerifyBamID](#VerifyBamID)
  * [Somalier](#Somalier)
  * [Merqury](#Merqury)
  * [sex-verify](#sex-verify)
* [Genome assembly](#Genome-assembly)
* [Genome alignment](#Genome-alignment)
* [Variant calling](#variant-calling)
* [SV merging](#SV-merging)
* [Annotation](#Annotattion)
* [Methylation](#methylation)
* [Housekeeping](HOUSEKEEPING.md)
* [Citation](#citation)

## Sample
### Sample origin/cohort
This study comprised 189 individuals (51 families) from the SSC, SAGE, and Rett-like cohorts, and the methodology is applicable to families with any rare disease.

| Count | Sex (proband-sibling) | Family type |
|:-----:|:---------------------:|------------:|
|  12   |          F-F          |        quad |
|  16   |          F-M          |        quad |
|   3   |          M-F          |        quad |
|   5   |          M-M          |        quad |
|  13   |           F           |        trio |
|   2   |           M           |        trio |

The sample manifest is available in the supplementary data of the publication.

[:arrow_double_up:](#table-of-contents)
## QC
##### back-reference-qc ([Kraken2](https://github.com/DerrickWood/kraken2))
* Use this pipeline to check for non-human contamination in reads.
  * Minimal requirement: FASTQ

##### [ntsm](https://github.com/JustinChu/ntsm)
* Use this tool/pipeline to assess inter-sample contamination.
    * Minimal requirement: FASTQ

##### [VerifyBamID](https://github.com/Griffan/VerifyBamID)
* Use this tool/pipeline to assess both non-human contamination and inter-sample contamination.
    * Minimal requirement: BAM

##### [Somalier](https://github.com/brentp/somalier)
* Use this tool/pipeline to assess inter-sample contamination, as well as ancestry and relatedness.
    * Minimal requirement: BAM

##### [Merqury](https://github.com/marbl/merqury)
* Use this tool/pipeline to assess genome assembly quality.
    * Minimal requirement: FASTQ and its own Illumina

##### [sex-verify](pipeline_scripts/sex-veriy)
* Use this to verify sex per cell or sample.
    * Minimal requirement: BAM
    * [click here for notes](notes/sex-verify.md)

[:arrow_double_up:](#table-of-contents)
## Genome assembly
This step produces FASTA files.
##### [hifiasm](https://github.com/chhylp123/hifiasm)
* Use this pipeline/tool to assemble sample genomes. Trio-phased assembly requires parental Illumina data as input.

* Version used for all our samples: hifiasm 0.16.1 with HiFi data only.

##### [fix-sex-chromosome](pipeline_scripts/fix-sex-chr)
* Use this pipeline to correct partially phased sex chromosomes in autism family fathers.

##### [Contiguous chromosome X/Y](https://github.com/projectoriented/contiguous-X)
* Use this pipeline to build contiguous sex chromosomes.

[:arrow_double_up:](#table-of-contents)
## Genome alignment
This step produces aligned BAM files. 
##### [pbmm2](https://github.com/PacificBiosciences/pbmm2)
* Use this pipeline to align HiFi FASTQ files to the reference genome. 

[:arrow_double_up:](#table-of-contents)
## Variant calling
##### [PAV](https://github.com/EichlerLab/pav)
* Use this tool to call SVs with assemblies. ([instructions](notes/pav.md))

##### [PBSV](https://github.com/PacificBiosciences/pbsv)
* Use this tool to call SVs with alignment (pbmm2 output).

##### [Sniffles](https://github.com/fritzsedlazeck/Sniffles)
* Use this tool to call SVs with alignment (pbmm2 output).

[:arrow_double_up:](#table-of-contents)
## SV merging
These steps are performed using [Truvari](https://github.com/ACEnglish/truvari) in sequential order.

#### 1. Intra-sample merge.

```shell
bcftools merge --thread {threads} --merge none --force-samples -O z -o {output.vcf.gz} {input.vcf1.gz} {input.vcf2.gz} {input.vcf3.gz}
truvari collapse -i {input.vcf.gz} -c {output.removed.vcf.gz} --sizemin 0 --sizemax 1000000 -k maxqual --gt het --intra --pctseq 0.90 --pctsize 0.90 --refdist 500 | bcftools sort --max-mem 8G -O z -o {output.collapsed.vcf.gz}
```

#### 2. Inter-sample merge.
```shell
bcftools merge --threads {threads} --merge none --force-samples --file-list {input.vcflist} -O z | bcftools norm --threads 15 --do-not-normalize --multiallelics -any --output-type z -o {output.mergevcf.gz}
truvari collapse --input {input.mergevcf.gz} --collapsed-output {output.removed_vcf.gz} --sizemin 0 --sizemax 1000000 --pctseq 0.90 --pctsize 0.90 --keep common --gt all | bcftools sort --max-mem {resources}G --output-type z > {output.collapsed_vcf.gz}
```

#### 3. [Rare SV pool discovery](pipeline_scripts/rareSVpool) of [an example input](https://eichlerlab.gs.washington.edu/public/rareSVpool/example_files).
```shell
python rareSVpool.py {input.collapsed_sv}
```
#### 4. De novo validation
* Initial caller support using [Truvari](https://github.com/ACEnglish/truvari)
* Callable region evaluation using [BoostSV](https://github.com/jiadong324/BoostSV)
* Genotyping support using [kanpig](https://github.com/ACEnglish/kanpig)
* Rare TR expansions/contractions using [TRGT](https://github.com/PacificBiosciences/trgt)
* Multiple sequence alignment (MSA) using [MAFFT](https://github.com/GSLBiotech/mafft)
* Read-based support validation using [subseq](https://github.com/EichlerLab/subseq) or [notes here](notes/denovo_correction.md)
* Manual inspection using IGV

[:arrow_double_up:](#table-of-contents)
## Annotation (GRCh38)
* Gene and location annotation using [AnnotSV](https://github.com/lgmgeo/AnnotSV), and then simplified by using [sim_annotSV.py](pipeline_scripts/comREG/sim_annotSV.py)
* CADD score using [CADD-SV](https://github.com/kircherlab/CADD-SV)
* Regulatory annotation using [REG data](https://eichlerlab.gs.washington.edu/public/comREG/data/) and [comREG](pipeline_scripts/comREG/regulation.snakefile)
* Combine all annotations from [comREG](pipeline_scripts/comREG/)

[:arrow_double_up:](#table-of-contents)
## [Methylation](https://github.com/projectoriented/continuous-methylation)
This step produces methylation bed files and corresponding bigwig files.

[:arrow_double_up:](#table-of-contents)
## Citation
For citation, please refer to our paper at: https://www.medrxiv.org/content/10.1101/2025.07.21.25331932v1

[:arrow_double_up:](#table-of-contents)

