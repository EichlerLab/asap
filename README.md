# ASAP
Autism Susceptibility Analysis Pipeline with a focus on Structural Variants (SVs). This repository documents the tasks involved in this project, which may be executed either sequentially or asynchronously. The approach used for rare variant or pathogenic candidate discovery in this study can be applied broadly to families (or individuals) affected by any rare disease.

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
The batch 1 of the study comprised 189 individuals (51 families) from the SSC, SAGE, and Rett-like cohorts, and the methodology is applicable to families (or individuals) with any rare disease.

The sample manifest is available in the supplementary data of the publication.

[:arrow_double_up:](#table-of-contents)
## QC
##### [back-reference-qc](https://github.com/EichlerLab/back-reference-qc)
* Use this pipeline (Kraken2) to check for non-human contamination in reads.
  * Minimal requirement: FASTQ

##### [sample-id-check](https://github.com/EichlerLab/sample-id-check)
* Use this pipeline (NTSM and VerifyBamID) to assess inter-sample contamination and verify sample identity.
    * Minimal requirement: FASTQ

##### [Somalier](https://github.com/brentp/somalier)
* Use this tool to assess inter-sample contamination, as well as ancestry and relatedness.
    * Minimal requirement: BAM

##### [assembly_qc](https://github.com/EichlerLab/assembly_qc)
* Use this pipeline to assess genome assembly quality.
    * Minimal requirement: FASTQ and its own Illumina

##### [sex-verify](pipeline_scripts/sex-veriy)
* Use this tool to verify the presence of the Y chromosome in each cell or sample.
    * Minimal requirement: BAM
    * [click here for notes](notes/sex-verify.md)

[:arrow_double_up:](#table-of-contents)
## Genome assembly
This step produces FASTA files.
##### [hifiasm](https://github.com/chhylp123/hifiasm)
* Use this tool to assemble sample genomes. Trio-phased assembly requires parental Illumina data as input.

* Version used across samples: hifiasm 0.16.1 or 0.25.0 with HiFi data only (samples were analyzed in stages).

##### [fix-sex-chromosome](pipeline_scripts/fix-sex-chr)
* Use this pipeline to correct partially phased sex chromosomes in autism family fathers, ensuring that hap1 corresponds to the Y chromosome and hap2 corresponds to the X chromosome.

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
python rareSVpoolv2.py --input {input.collapsed_sv} --sample {input.sample_manifest} --outdir {output.rare} {optional: --freq 1}
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
## Annotation
* Gene and location annotation (GRCh38) using [AnnotSV](https://github.com/lgmgeo/AnnotSV), and then simplified by using [sim_annotSV.py](pipeline_scripts/comREG/sim_annotSV.py) .
* CADD score (GRCh38) using [CADD-SV](https://github.com/kircherlab/CADD-SV) .
* Combined Regulatory Annotation (GRCh38/CHM13) using [REG data](https://eichlerlab.gs.washington.edu/public/comREG/data/) and [comREG](pipeline_scripts/comREG/) .


[:arrow_double_up:](#table-of-contents)
## [Methylation](https://github.com/projectoriented/continuous-methylation)
This step produces methylation bed files and corresponding bigwig files.

[:arrow_double_up:](#table-of-contents)
## Citation
For citation, please refer to our paper at: https://www.nature.com/articles/s41467-026-68378-4
[:arrow_double_up:](#table-of-contents)

## Contact
Questions, feedback, or bug reports are always welcome — reach out to Yang Sui at yangsui@uw.edu.
