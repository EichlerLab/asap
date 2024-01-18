# SV-Pop
https://github.com/EichlerLab/svpop
## Directory setup
```commandline
ln -s /net/eichler/vol28/software/pipelines/svpop/svpop-3.4.2/runsnake
ln -s /net/eichler/vol28/software/pipelines/svpop/svpop-3.4.2/runlocal
```

### Contents of config/config.json
```json
{
    "variant_table": "config/samples.tab",
    "reference": "/net/eichler/vol28/eee_shared/assemblies/hg38/no_alt/hg38.no_alt.fa",
    "ucsc_ref_name": "hg38",
    "samplelist": {
        "family": [
            "11071_fa",
            "11071_mo",
            "11071_p1",
            "11071_s1",
            "11080_fa",
            "11080_mo"
        ]
    },
    "callerset": {},
    "sampleset": {
        "asd": {
            "sourcetype": "caller",
            "sourcename": "pav-hifi",
            "merge": {
                "svindel:ins,del,insdel": "nr::szro(0.5,200,4)"
            },
            "name": "ASD quads",
            "description": "Reflect aggregated and filtered PAV + pbsv variants onto the healthy HPRC + 1000g populations"
        }
    },
    "merge_def": {
        "szro-50-200": "nr::szro(0.5,200,4)",
        "szro-80": "nr::szro(0.8,,4):match(0.8)",
        "sv-exact": "nr::exact"
    }
}
```

### Contents of config/samples.tab
```shell
NAME    TYPE    SAMPLE  DATA    VERSION
pbsv-hifi       pbsv    DEFAULT GRCh38/variant_call/pbsv/{sample}/pbsv_{sample}_{vartype}.vcf.gz     2.9.0
pav-hifi        pavbed  DEFAULT GRCh38/pav/1.1.2/all/results/{sample}/bed/ 1.1.2
```

## Run
```shell
# getting the bed files for each sample
./runsnake 80 results/variant/caller/{pav-hifi,pbsv-hifi}/${sample}/{all,lc}/{all,notr}/bed/sv_{ins,del,insdel}.bed.gz

# Merge the variants in sampleset
echo results/variant/sampleset/asd/family/{all,lc}/{all,notr}/bed/sv_{insdel,ins,del}.bed.gz | tr ' ' '\n'
```