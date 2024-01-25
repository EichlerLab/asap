# Verify the sex of your sample

Get the resources [here](../pipeline_scripts/sex-verify).
## Directory setup
```txt
.
├── config.yaml
├── manifest.tab
├── qc-sex.smk
├── run.sh
└── sex-estimator.py
```

## Contents of manifest.tab
```tsv
sample	file_path
11080_fa	/net/eichler/vol28/projects/long_read_archive/nobackups/clinical/11080_fa/raw_data/PacBio_HiFi/fofn/ccs/fastq.fofn
```

If you want to check for ONT data, then just make a manifest called manifest-ont.tab or whatever you'd like. You just need to let the pipeline know.
## Run
```shell
# example for ont
./work.sh 80 --config manifest=manifest-ont.tab tech=ont -p

# example for ont
./work.sh 80 -p
```