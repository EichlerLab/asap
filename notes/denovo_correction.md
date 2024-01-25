# De novo correction
This document walks through the correction of discovery curve outputs regarded as singletons.


### Set up the directory & modify definitions.snakefile file
```shell
cd /net/eichler/vol28/projects/autism_genome_assembly/nobackups/post_processing/GRCh38/denovo_correction

mkdir -p beds

module load miniconda/4.12.0
cd beds && ./get_beds.py asd-subset.tsv.gz

# make sure to modify definitions.snakefile if necessary
less definitions.snakefile
``` 

### Get the targets & run
```shell
ls -d beds/* | grep -Eo "[BK0-9]{5}_.{2}" | sort -u | while read line; do ./get_targets.sh $line; done > targets.txt
./runsnake 100 $(cat targets.txt) -p
```

### Get the corrected output
```shell
./correct_denovo.sh
less -S denovo_corrected.tsv.gz
```

### Merge with original
Please merge the output of 1) data table (the filtered one if it exists), 2) `asd-subset.tsv.gz`, and 3) `denovo_corrected.tsv.gz` based on the IDs. This is a good opportunity to customize the final SV table you view.