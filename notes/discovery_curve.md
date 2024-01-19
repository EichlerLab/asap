## Directory setup
```commandline
ln -s /net/eichler//vol28/software/pipelines/disc_curve/runsnake
ln -s /net/eichler//vol28/software/pipelines/disc_curve/runlocal
```

### Contents of config.yaml
```yaml
# the starting backgrounds
background: /net/eichler/vol28/projects/hprc/nobackups/data_table/qc/hprc_hgsvc/tsv/variants_hprc_hgsvc_sv_insdel.tsv.gz
sample_order: ["14694_s1","14616_p1"]
full_curve: /net/eichler/vol28/projects/hprc/nobackups/data_table/qc/hprc_hgsvc/tables/class_count/sv/class_count_all_sv_insdel_all.tsv.gz
```

### Contents of manifest.tab
```tsv
SAMPLE DATA_TABLE PAV_DIR MERGE_MAP
14694_s1 /net/eichler/vol28/projects/autism_genome_assembly/nobackups/post_processing/GRCh38/data_table/tsv/variants_asd_families_sv_insdel_filt.tsv.gz /net/eichler/vol28/projects/autism_genome_assembly/nobackups/variant_calling/GRCh38/pav/1.1.2/all/ /net/eichler/vol27/projects/autism_genome_assembly/nobackups/post_processing/GRCh38/data_table/sections/asd_families/base_table/merge_map_sv_insdel.tsv.gz
```

## Run
```shell
./runlocal 30 -p
```

### When it is done
Use the script [`process_DC_output.py`](../pipeline_scripts/process_DC_output.py)
```shell
./process_DC_output.py \
  --merge_map data_table/sections/asd_families/base_table/merge_map_sv_insdel.tsv.gz \
  --disc_curve discovery_curve/results/variants_bg-ALL-sv_insdel.tsv.gz \
  -t 20 \
  -o all.tsv.gz
```