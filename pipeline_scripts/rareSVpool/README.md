# rareSVpool

Genome-wide rare structural variant (SV) pool discovery for LRS-based rare disease studies. Given a cohort of collapsed/merged SV calls (e.g. Truvari `collapse` output) and a sample manifest describing family structure, the script identifies SVs that are specific to children/case and absent from unrelated controls, with het/hom-aware and sex-aware handling of chrX/chrY — filtering out sex-linked variants that also occur in same-sex control individuals.


## Requirements

- Python 3
- `pandas`, `numpy`

## Usage

```bash
python rareSVpoolv2.py \
  --input {input.disco_truvari_collapsed.txt.gz} \
  --sample {input.samples.txt} \
  --outdir {output.rare} \
  {optional: --freq 1}
```

### Arguments

| Flag | Required | Description |
|---|---|---|
| `--input` | yes | Collapsed SV table (tab-separated), e.g. Truvari collapse output. Must include `#CHROM`, `ID`, `MERGE_SAMPLES`, `MERGE_GT`. |
| `--sample` | yes | Sample sheet (tab-separated) with columns `sample`, `sex`, `group`, `famid`. Unknown sex can be estimated with `sex-estimator.py`. |
| `--outdir` | yes | Output directory. |
| `--freq` | no | Keep rare het/hom SVs present in at most this many families. If omitted, no frequency filtering is applied — **recommended** for comprehensive pathogenic SV discovery, so no candidate is excluded on frequency alone.|

### Sample sheet (`--sample`) format

Tab-separated with one row per sample:

- `sample` — sample ID including both controls and query samples, matching the IDs used in `MERGE_SAMPLES`.
- `sex` — `M` or `F`.
- `group` — normalized to one of `control`, `father`, `mother`, `proband`, `sibling`. Accepted aliases:
  - control: `control`, `ctrl`, `ctr`
  - father: `father`, `fa`, `dad`
  - mother: `mother`, `mo`, `mom`
  - proband: `proband`, `pro`, `p1`, `p2`
  - sibling: `sibling`, `sib`, `s1`, `s2`
- `famid` — family ID, normalized to the proband(case)'s ID to recognize potential duo/trio/quad family structures.

### Input SV table format

Tab-separated, one row per SV, expected to include at minimum:

- `#CHROM`, `ID`
- `MERGE_SAMPLES` — comma-joined list of sample IDs carrying the SV.
- `MERGE_GT` — comma-joined list of genotypes (e.g. `0|1`, `1|1`), aligned positionally with `MERGE_SAMPLES`.

## Output files

Written to `--outdir`:

| File | Contents |
|---|---|
| `disco_truvari_collapsed_processed.txt.gz` | Full input table with added sample/sex/class annotation columns, after removing mis-called chrY rows. |
| `disco_truvari_collapsed_processed_rareSVpool.txt` | The rare SV pool, expanded to one row per carrying child, with `Category`, `FAMID`, `role`, `Total_N_Fam`, `Total_N_Fam_homo`. |
| `disco_truvari_collapsed_processed_rareSVpool_ct.txt` | Count of rare SVs per `Sample`/`role`/`Sex`. |
| `disco_truvari_collapsed_processed_rareSVpool{freq}fam.txt` | (only if `--freq` set) Pool subset where family occurrence ≤ `freq`. |
| `disco_truvari_collapsed_processed_rareSVpool{freq}fam_ct.txt` | (only if `--freq` set) Per-sample counts for the filtered subset. |

The script also prints to stdout: total/original SV counts, mean rare SVs per child, and (if `--freq` is set) the mean rare SVs and mean singletons per child at that frequency threshold.

## Notes

- Male genotypes on chrX are corrected from `1|1` to `0|1` (hemizygous) before rarity is assessed.
- A variant is only called "rare" on chrX if no *sex-matched* control carries it — a male proband's chrX call is not disqualified by a female control carrying the same call, and vice versa.
- Two rounds of chrY filtering guard against dirty merges: an initial pass drops chrY SVs with no male carrier at all, and a later pass drops chrY calls contaminated by female carriers or mixed-sex controls, resulting in `total/original SV counts`.
- `Total_N_Fam` / `Total_N_Fam_homo` exclude female samples from chrY family counts, since females cannot carry Y-chromosome variants.
- `--freq` is generally **not recommended** for complete pathogenic SV discovery. Simple frequency filtering doesn't account for ONT/low-coverage SV calling quality (genotypes and breakpoints are noisier than HiFi calls) or for sex-specific zygosity on chrX — e.g. a hemizygous call in a male proband and a heterozygous call in a female control represent different biological states but would both simply count toward "occurrence," conflating carriers that the sex-matched filtering logic above is specifically designed to keep apart. Leaving `--freq` unset keeps the full, sex-aware pool intact for review and downstream verification.
- This script is under active development; logic, defaults, and output formats may change. This version is similar to the one used in the autism study (https://www.nature.com/articles/s41467-026-68378-4), generalized for broader use across family structures and study designs. 

## Contact

Questions, feedback, or bug reports are always welcome — reach out to Yang Sui at yangsui@uw.edu.
