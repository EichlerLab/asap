# Combined Regulatory Annotation (comREG)

Snakemake pipeline that annotates a set of SVs/regions against a panel of regulatory, functional, and repeat-element tracks on **GRCh38**. A **CHM13** version of this pipeline (`regulation_CHM13.snakefile`) is also available, using the equivalent CHM13-based annotation tracks with the same logic.

## What it does

1. Intersects the input variant BED file against every annotation track listed in `manifest.txt` using `bedtools intersect`, collapsing multiple overlapping features per variant into a single semicolon-joined value.
2. Runs a few tracks with track-specific logic: `platinumTRs` (tandem repeats) requires the variant to be overlapped by at least 50% of its length (`-f 0.5`); non-coding constraint and segmental-duplication tracks keep only the single highest-scoring overlap per variant.
3. Merges every per-track result into one wide table indexed by variant `ID`.
4. Parses the GENCODE annotation string into gene lists by region (`CDS`, `5UTR`, `3UTR`, and all genes overall).
5. Derives summary flags for whether a variant falls in any regulatory element (`REG`), any regulatory element plus TF clusters (`REGwTF`), non-brain regulatory elements only (`otherREG`), or brain-specific marks only (`brainREG`).
6. Assigns each variant a single `Location_GENCODE` call (`CDS` > `5UTR` > `3UTR` > `REG` > `.`) summarizing its most relevant genic/regulatory context.
7. Writes the final merged, annotated table to `allreg_sv.txt`.

## Requirements

- Snakemake
- `bedtools`
- Python 3 with `pandas`, `numpy`
- `miniconda/4.12.0` environment module (cluster-specific; adjust `envmodules` if running elsewhere)

## Inputs

| File | Description |
|---|---|
| `annotation_reg_input.bed` | Variant/SV BED file (tab-separated, with header) to annotate. Must include an `ID` column (4th column, used as the join key throughout) and coordinate columns compatible with `bedtools intersect`. |
| `manifest.txt` | Tab-separated manifest with columns `peak` (track name, used as index) and `path` (file path to that track's BED file). Drives the generic per-track intersection rule. Manifest track files are available for download at https://eichlerlab.gs.washington.edu/public/comREG/data/. |


## Output

`allreg_sv.txt` — one row per variant `ID`, with:

- Per-track annotation columns (gene/feature names overlapping the variant, `;`-joined if multiple).
- `gene`, `CDS_genes`, `5UTR_genes`, `3UTR_genes` — parsed from the `GENCODE` column.
- `REG`, `REGwTF`, `otherREG`, `brainREG` — `YES`/`NO` summary flags described above.
- `Location_GENCODE` — single-label genic/regulatory context call per variant.

## Notes

- Gene-level constraint and neurodevelopmental/ASD gene-list joins (pLI, LOEUF, SFARI, Fu NDD664, Wang NDD615, Satterstrom ASD102, Sui NDD810, MRG, CMRG, HPO, phenotype) are present in the code but **commented out**, since they point to lab-internal file paths. Uncomment and repoint them to your own copies of these resources to re-enable that layer of annotation.
- The input BED filename (`annotation_reg_input.bed`), manifest filename (`manifest.txt`), and the header/column names the script expects (e.g. `ID`, `GENCODE`) are **hard-coded** in the Snakefile rather than configurable. Match these names exactly, or edit the script if your files use different names.
- This pipeline is under active development; tracks, defaults, and output columns may change.

## Contact

Questions, feedback, or bug reports are always welcome — reach out to Yang Sui at yangsui@uw.edu.
