# ASAP
Autism Susceptibility Analysis Pipeline with focus in Structural Variants.

This repo logs tasks involved for this project- either executed sequentially or asynchronously.  

##### Table of Contents
* [Inputs](#inputs)
* [QC](#qc)
  * [back-reference-qc](#back-reference-qc)
  * [merqury](#merqury)
  * [verifybamid](#verifybamid)
  * [sex-verify](#sex-verify)
* [Genome Assembly](#genome-assembly)
* [Alignment](#alignment)
* [Variant calling](#variant-calling)
* [Getting a refined call set](#refined-callset)
* [Methylation](#methylation)
* [Analyses](#analyses)
  * [Contiguous chromosome X](#contiguous-chromosome-x)
  * [DNM visualization](#dnm-visualization)
  * [DMR identification](#dmr-identification)
* [FAQ](#faq)
* [Paths](#paths)
* [Housekeeping](HOUSEKEEPING.md)

## Inputs
### Sample origin/cohort
This pertains to SSC + SAGE cohort.

| Count | Sex (proband-sibling) | Family type |
|:-----:|:---------------------:|------------:|
|  11   |          F-F          |        quad |
|  17   |          F-M          |        quad |
|   4   |          M-F          |        quad |
|   5   |          M-M          |        quad |
|  10   |           F           |        trio |
|   2   |           M           |        trio |
|   1   |         F-F-F         |       penta |

Sample sheet is here: `/net/eichler/vol28/projects/autism_genome_assembly/nobackups/sample_info.tab`

[:arrow_double_up:](#table-of-contents)
## QC
##### back-reference-qc
* use this pipeline for non-human contamination of reads.
  * minimal requirement: fastq.gz
  * https://eichlerlab.gs.washington.edu/help/wiki/doku.php?id=users:lettucerap:back_reference_qc
##### merqury
* use this tool/pipeline to assess quality of genome assembly
    * minimal requirement: fastq.gz and its own Illumina
    * https://eichlerlab.gs.washington.edu/help/wiki/doku.php?id=users:merqury
##### verifybamid 
use this tool/pipeline to assess contamination of non-humanness as well as inter-sample contamination.
    * minimal requirement: bam
    * https://eichlerlab.gs.washington.edu/help/wiki/doku.php?id=users:nidhi12k:vbi
##### sex-verify
use this to check either per cell or sample for sex verification
    * minimal requirement: bam
    ```shell
    cd /net/eichler/vol28/projects/autism_genome_assembly/nobackups/qc/sex-verify
    # modify the config.yaml and manifest.tab as needed
    ./run.sh 30 -p
    ```

[:arrow_double_up:](#table-of-contents)
## Genome Assembly
This step produces a fasta file.
* hifiasm: use this pipeline/tool to assemble sample genome
  * **trio-phased requires** parental Illumina data as input
  ```shell
  cd /net/eichler/vol28/projects/autism_genome_assembly/nobackups/assemblies/hifiasm_0.16.1 && ./run_hifiasm.sh 30 -p
  ```
* Version used for all our samples ATM (Jan 16, 2024): hifiasm 0.16.1 with just HiFi data.

[:arrow_double_up:](#table-of-contents)
## Alignment
This step is produces a BAM. And can be achieved via: https://eichlerlab.gs.washington.edu/help/wiki/doku.php?id=users:wharvey:align_all for all.
* PacBio_HiFi fastq.gz input, the pipeline uses [pbmm2](https://github.com/PacificBiosciences/pbmm2)
* Illumina fastq.gz, the pipeline uses bwa2 and you can find fofn here:
    ```shell
    ls -lrtha /net/eichler/vol28/projects/autism_genome_assembly/nobackups/data/Illumina/WGS/fastq/fofn
    ```

[:arrow_double_up:](#table-of-contents)

## Variant calling
* [pbsv](https://github.com/PacificBiosciences/pbsv): use pbmm2 output as input for this
  ```shell
  cd /net/eichler/vol28/projects/autism_genome_assembly/nobackups/variant_calling/GRCh38/pbsv && ./run_pbsv.sh 80 all
  # make sure you populate pbsv.tab
  ```
* [PAV](https://github.com/EichlerLab/pav): use hifiasm assembly output for this
  * [instructions](notes/pav.md)
  ```shell
  cd /net/eichler/vol28/projects/autism_genome_assembly/nobackups/variant_calling/GRCh38/pav/1.1.2 && ./run_pav.sh 80 all assemblies-all.tsv
  # make sure you populate assemblies-all.tab
  ```

[:arrow_double_up:](#table-of-contents)

## Refined callset
The steps here are by sequential order.

#### SUBSEQ for SV validation
```shell
cd /net/eichler/vol28/projects/autism_genome_assembly/nobackups/post_processing/GRCh38/sub_seq
./get_targets.sh ${sample} > targets.txt
./runsnake 30 $(cat targets.txt)
```
#### SVPOP sampleset merging of SVs
* [SVPOP](https://github.com/EichlerLab/svpop): [instructions here](notes/svpop.md)
1. Get the reformatted bed files for each caller, e.g. pav-hifi
2. Get the sampleset merge bed file
3. Get the annotations
```shell
cd /net/eichler/vol28/projects/autism_genome_assembly/nobackups/post_processing/GRCh38/svpop
./runsnake 100 -p $(./get_intersect_targets.sh 14455_p1)
./runsnake 100 $(cat ./get_anno_targets.sh asd family | grep -v "inv") -p
```

#### Data table to summarize SVPOP annotations
```shell
cd /net/eichler/vol28/projects/autism_genome_assembly/nobackups/post_processing/GRCh38/data_table
./runsnake 30 -p $(cat first_target.txt) && ./runsnake 30 -p $(cat second_target.txt)
.filter_data_tbl.py tsv/variants_asd_families_sv_insdel.tsv.gz
./get_supported_variants.py
```

#### Discovery curve
Instructions are [here](notes/discovery_curve.md)

Highlights:
* `config.yaml`:
  * `sample_order`: desired order for the output cumulative figure
  * `full_curve`: hgsvc + hprc cumulative counts of the truthset
  * `background`: truth set data table of SV insdel.
* `manifest.tab`:
  * `DATA_TABLE`: data table output that contains samples within the `sample_order`
  * `PAV_DIR`: run directory of pav, the pipeline searches for this file `results/${sample}/bed/sv_{ins,del}.bed.gz` from the run directory- and the critical columns it searches for are ID,#CHROM,POS,END,SVTYPE,SVLEN
  * `MERGE_MAP`: this can be found in the output of data table. example output here: `/net/eichler/vol28/projects/autism_genome_assembly/nobackups/post_processing/GRCh38/data_table/sections/asd_families/base_table/merge_map_sv_insdel.tsv.gz`

```shell
cd /net/eichler/vol28/projects/autism_genome_assembly/nobackups/post_processing/GRCh38/discovery_curve
# modify manifest.tab and config.yaml
./runlocal -p
```
#### DNM validation
Instructions are [here](notes/denovo_correction.md)

[:arrow_double_up:](#table-of-contents)

## Methylation
This step produces a methylation bed file and bigwig files of the beds. The pipeline for both technologies are here: https://eichlerlab.gs.washington.edu/help/wiki/doku.php?id=users:lettucerap:methylation
* When input is HiFi data, [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools) is used.
* When input is ONT data, [modkit](https://github.com/nanoporetech/modkit) is used.

[:arrow_double_up:](#table-of-contents)

## Analyses
#### Contiguous chromosome X
* https://github.com/projectoriented/contiguous-X

#### DNM visualization
* https://eichlerlab.gs.washington.edu/help/wiki/doku.php?id=users:lettucerap:ava_svbyeye

#### DMR identification
Group-comparison pipeline
* https://eichlerlab.gs.washington.edu/help/wiki/doku.php?id=users:lettucerap:dss_snakemake

## FAQ
1. I think the sample is contaminated, how should I investigate?
   1. Try verifybamid either for each cell *or* for the entire sample aligned against the same reference as its Illumina bam.
   2. Try back-reference-qc to see how much non-human contamination there is and remove it from the fastq.gz
2. I think the sample is of different sex, what to do?
   1. Try looking at coverage over the SRY gene for quick visualization
   2. Estimate the sex using [this script](https://github.com/projectoriented/bio-utils/blob/main/sex-estimator.py) either for each cell or the entire sample. The input must be a BAM. FYI: the script is already in one of the QC folders but it may not be the newest version.
3. How do I make ONT fofns?
   1. https://eichlerlab.gs.washington.edu/help/wiki/doku.php?id=users:lettucerap:make_ont_fofn
4. Where are the HiFi fofns?
   1. $LRA/clinical/${sample}/raw_data/PacBio_HiFi/fofn/ccs/fastq.fofn

[:arrow_double_up:](#table-of-contents)

## Paths
```shell
LRA=/net/eichler/vol28/projects/long_read_archive/nobackups
sample=14455_p1

# File of File Names (fofn)
# hifi
ls -lrtha $LRA/clinical/${sample}/raw_data/PacBio_HiFi/fofn/ccs/fastq.fofn

# ont- get link from FAQ for instructions to generate independently
cd /net/eichler/vol28/projects/autism_genome_assembly/nobackups/methylation/fastq_fofn/ && ./get-fastq.sh

cd /net/eichler/vol28/projects/autism_genome_assembly/nobackups/methylation/unmapped_bam_fofn/ && ./get-unmapped-bam.sh

# hifiasm assemblies
ls -lrtha $LRA/clinical/${sample}/assemblies/hifiasm
```

[:arrow_double_up:](#table-of-contents)
