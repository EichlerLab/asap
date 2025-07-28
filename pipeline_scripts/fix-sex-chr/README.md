# reassign-sex-chr

This pipeline was first designed by @projectoriented and modified by Yang Sui (yangsui@uw.edu)

## Methodology

All contigs which align to the sex chromosomes are assigned into hap1 (chrY) and hap2 (chrX) based on alignments to a reference

1. Alignments are made to a given reference with chrX and chrY (in most cases CHM13)
2. Contigs are sorted into one of two categories:
    - Complement (Split)
      - Those contigs which do not have an overlap with another contig in the opposite haplotype but are located in the haplotype not designated to that sex chromosome
    - Align (Duplicated)
      - Those contigs which have overlapping alignments with contigs in another haplotype for the same chromosome and coordinates
    -  By default these criteria are required to reach a 95% reciprocal overlap to determine similarity
3. The pipeline can either produce just a list of contigs and their desired assignments, or reorder the input fasta files accordingly

## Running the pipeline
### Step 1. Prepare directory
```
config
├── config.yaml
└── manifest.tab
runsnake
```
### Step 2. Prepare inputs
#### config/config.yaml
```
reference: T2T-CHM13v2.fasta

#### OPTIONAL ARGS
minimap_params: '-x asm20 --secondary=no -s 25000'
manifest: config/manifest.tab
```
#### config/manifest.tab
```bash
sample  hap1    hap2
your_sample path/to/hap1_asm    /path/to/hap2_asm
```
#### Step 3. Start the analysis!
```bash
./runsnake 30
```

