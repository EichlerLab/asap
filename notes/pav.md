# PAV

## Directory setup
```commandline
ln -s /net/eichler/vol28/projects/structural_variation/nobackups/pipelines/pav/2.3.4/runsnake
ln -s /net/eichler/vol28/projects/structural_variation/nobackups/pipelines/pav/2.3.4/runlocal
```

### Contents of config.json
```json
{
    "reference": "/net/eichler/vol28/eee_shared/assemblies/hg38/no_alt/hg38.no_alt.fa",
    "assembly_table": "manifest.tab"
}
```

### Contents of manifest.tab (a file of sample paths and names)
```commandline
NAME    HAP1    HAP2
11918_s1        11918_s1.hifiasm.dip.hap1.p_ctg.gfa.fasta      11918_s1.hifiasm.dip.hap2.p_ctg.gfa.fasta
11918_p1        11918_p1.hifiasm.dip.hap1.p_ctg.gfa.fasta      11918_p1.hifiasm.dip.hap2.p_ctg.gfa.fasta
```

## Run
```shell
module load miniconda/4.12.0
./runsnake 30 -np # this is dry run
./runsnake 30 -p # actually submitting jobs
```
