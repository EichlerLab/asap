# De novo correction
This document walks through the correction of discovery curve outputs regarded as singletons.


### Set up the directory & modify definitions.snakefile file
```shell
cd /net/eichler/vol28/projects/autism_genome_assembly/nobackups/post_processing/GRCh38/denovo_correction

mkdir -p beds

module load miniconda/4.12.0
cd beds && ./get_beds.py file-after-discovery-curve-and-appended-family-info.tsv.gz

# make sure to modify definitions.snakefile if necessary
less definitions.snakefile
``` 

### Get the targets & run
```shell
ls -d beds/* | grep -Eo "[BK0-9]{5}_.{2}" | sort -u | while read line; do ./get_target.sh $line; done > targets.txt
./runsnake 100 $(cat targets.txt) -p
```

### Get the corrected output
```shell
./correct_denovo.sh
less -S denovo_corrected.tsv.gz
```

### Merge with original
```shell

```