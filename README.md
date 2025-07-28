# ASAP
Autism Susceptibility Analysis Pipeline with focus in Structural Variants.

This repo logs tasks involved for this project- either executed sequentially or asynchronously.  

##### Table of Contents
* [Inputs](#inputs)
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

## Inputs
### Sample origin/cohort
This pertains to SSC + SAGE + Rett-like cohort (189 individuals from 51 families).

| Count | Sex (proband-sibling) | Family type |
|:-----:|:---------------------:|------------:|
|  12   |          F-F          |        quad |
|  16   |          F-M          |        quad |
|   3   |          M-F          |        quad |
|   5   |          M-M          |        quad |
|  13   |           F           |        trio |
|   2   |           M           |        trio |

Sample sheet is here: `/net/eichler/vol28/projects/autism_genome_assembly/nobackups/sample_info.tab`

`Dataset S1`

[:arrow_double_up:](#table-of-contents)
## QC
##### back-reference-qc ([Kraken2](https://github.com/DerrickWood/kraken2))
* use this pipeline for non-human contamination of reads.
  * minimal requirement: fastq.gz
  * [Internal path](https://eichlerlab.gs.washington.edu/help/wiki/doku.php?id=users:lettucerap:back_reference_qc)

##### [ntsm](https://github.com/JustinChu/ntsm)
* use this tool/pipeline to assess inter-sample contamination.
    * minimal requirement: fastq.gz

##### [VerifyBamID](https://github.com/Griffan/VerifyBamID)
* use this tool/pipeline to assess contamination of non-humanness as well as inter-sample contamination.
    * minimal requirement: bam
    * [Internal path](https://eichlerlab.gs.washington.edu/help/wiki/doku.php?id=users:nidhi12k:vbi)

##### [Somalier](https://github.com/brentp/somalier)
* use this tool/pipeline to assess inter-sample contamination as well as ancestry and relatedness.
    * minimal requirement: bam

##### [Merqury](https://github.com/marbl/merqury)
* use this tool/pipeline to assess quality of genome assembly
    * minimal requirement: fastq.gz and its own Illumina
    * [Internal path](https://eichlerlab.gs.washington.edu/help/wiki/doku.php?id=users:merqury)

##### [sex-verify](pipeline_scripts/sex-veriy)
* use this to check either per cell or sample for sex verification
    * minimal requirement: bam
    * [click here for notes](notes/sex-verify.md)

[:arrow_double_up:](#table-of-contents)
## Genome assembly
This step produces a fasta file.
* [hifiasm](https://github.com/chhylp123/hifiasm): use this pipeline/tool to assemble sample genome (trio-phased requires parental Illumina data as input).

* Version used for all our samples: hifiasm 0.16.1 with just HiFi data.

* [fix-sex-chromosome](pipeline_scripts/fix-sex-chr): use this pipeline to fix partially phased autism family fathers.
  * [Internal path](https://eichlerlab.gs.washington.edu/help/wiki/doku.php?id=users:lettucerap:sex_chromosome_grouping)

* [Contiguous chromosome X/Y](https://github.com/projectoriented/contiguous-X)
  * [Internal path](https://eichlerlab.gs.washington.edu/help/wiki/doku.php?id=users:lettucerap:contiguous-x)


[:arrow_double_up:](#table-of-contents)
## Genome alignment
This step is produces a BAM. And can be achieved via: [Internal path](https://eichlerlab.gs.washington.edu/help/wiki/doku.php?id=users:wharvey:align_all) 
* PacBio_HiFi fastq.gz input, the pipeline uses [pbmm2](https://github.com/PacificBiosciences/pbmm2)


[:arrow_double_up:](#table-of-contents)

## Variant calling
* [PAV](https://github.com/EichlerLab/pav): use hifiasm assembly output for this.([instructions](notes/pav.md))
* [pbsv](https://github.com/PacificBiosciences/pbsv): use pbmm2 output as input for this.
* [Sniffles](https://github.com/fritzsedlazeck/Sniffles): use pbmm2 output as input for this.


[:arrow_double_up:](#table-of-contents)

## SV merging
The steps here are using [Truvari](https://github.com/ACEnglish/truvari) by sequential order.

#### 1. sampleset merge.

```shell
bcftools merge --thread {threads} --merge none --force-samples -O z -o {output.vcf.gz} {input.vcf1.gz} {input.vcf2.gz} {input.vcf3.gz}
truvari collapse -i {input.vcf.gz} -c {output.removed.vcf.gz} --sizemin 0 --sizemax 1000000 -k maxqual --gt het --intra --pctseq 0.90 --pctsize 0.90 --refdist 500 | bcftools sort --max-mem 8G -O z -o {output.collapsed.vcf.gz}
```

#### 2. inter-sample nmerge.
```shell
bcftools merge --threads {threads} --merge none --force-samples --file-list {input.vcflist} -O z | bcftools norm --threads 15 --do-not-normalize --multiallelics -any --output-type z -o {output.mergevcf.gz}
truvari collapse --input {input.mergevcf.gz} --collapsed-output {output.removed_vcf.gz} --sizemin 0 --sizemax 1000000 --pctseq 0.90 --pctsize 0.90 --keep common --gt all | bcftools sort --max-mem {resources}G --output-type z > {output.collapsed_vcf.gz}
```

#### 3. [Rare SV pool discovery](pipeline_scripts/rareSVpool).
```shell
python rareSVpool.py {input.collapsed_sv}
```
#### 4. de novo validation
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
* Regulatory annotation using [comREG](pipeline_scripts/comREG/regulation.snakefile)
* Combine all annotations using [comREG](pipeline_scripts/comREG/combine_anno.py)

[:arrow_double_up:](#table-of-contents)

## [Methylation](https://github.com/projectoriented/continuous-methylation)
This step produces a methylation bed file and bigwig files of the beds. [Internal path](https://eichlerlab.gs.washington.edu/help/wiki/doku.php?id=users:lettucerap:methylation)


[:arrow_double_up:](#table-of-contents)


