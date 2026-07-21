#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import os


def args():
    p = argparse.ArgumentParser(description="Genome-wide rare SV pool discovery")
    p.add_argument("--input", required=True, help="Input collapsed SVs")
    p.add_argument("--sample", required=True, help="Sex and group information of the samples (unknown sex can be estimated using sex-estimator.py)")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--freq", required=False, type=int, default=None, help="Filter by occurrence (family count): keep rare het/hom SVs present in up to this many families. Optional.")
    return p.parse_args()


# ---- group vocabulary --------------------------------------------------------
# The --sample file must have columns: sample, sex, group, famid
# `group` is normalized to one of: control / father / mother / proband / sibling
# `famid` is normalized to the proband
# `sex` is important for sex-matched filtering
CONTROL_ALIASES = {"control", "ctrl", "ctr"}
FATHER_ALIASES  = {"father", "fa", "dad"}
MOTHER_ALIASES  = {"mother", "mo", "mom"}
PROBAND_ALIASES = {"proband", "pro", "p1", "p2"}
SIBLING_ALIASES = {"sibling", "sib", "s1", "s2"}

CHILD_GROUPS  = {"proband", "sibling"}
PARENT_GROUPS = {"father", "mother"}
FAMILY_GROUPS = CHILD_GROUPS | PARENT_GROUPS   # everything that is not a control


def canon_group(g):
    g = str(g).strip().lower()
    if g in CONTROL_ALIASES: return "control"
    if g in FATHER_ALIASES:  return "father"
    if g in MOTHER_ALIASES:  return "mother"
    if g in PROBAND_ALIASES: return "proband"
    if g in SIBLING_ALIASES: return "sibling"
    return g


def main():
    a = args()
    os.makedirs(a.outdir, exist_ok=True)

    ################ load sample metadata.
    df = pd.read_csv(a.input, sep='\t')
    meta = pd.read_csv(a.sample, sep='\t')

    sex_dict   = dict(zip(meta["sample"], meta["sex"]))
    group_dict = {s: canon_group(g) for s, g in zip(meta["sample"], meta["group"])}
    famid_dict = dict(zip(meta["sample"], meta["famid"]))

    def group_of(s):
        return group_dict.get(s, "")

    def famid_of(s):
        return famid_dict.get(s, "")

    # helpers over a comma-joined MERGE_SAMPLES / MERGE_1|1_All string
    def has_control(samples):
        if not isinstance(samples, str) or samples == "":
            return False
        return any(group_of(s) == "control" for s in samples.split(","))

    def has_child(samples):
        if not isinstance(samples, str) or samples == "":
            return False
        return any(group_of(s) in CHILD_GROUPS for s in samples.split(","))

    def has_family(samples):  # any non-control (parent or child) sample
        if not isinstance(samples, str) or samples == "":
            return False
        return any(group_of(s) in FAMILY_GROUPS for s in samples.split(","))

    ################ add sex info.
    def extract_samples_info(samples):
        sample_list = samples.split(",")
        all_sex = [sex_dict.get(s, "") for s in sample_list]
        ctr_samples    = [s for s in sample_list if group_of(s) == "control"]
        sample_samples = [s for s in sample_list if group_of(s) in CHILD_GROUPS]
        ctr_sex    = [sex_dict.get(s, "") for s in ctr_samples]
        sample_sex = [sex_dict.get(s, "") for s in sample_samples]
        return (
            ",".join(all_sex) if all_sex else "",
            ",".join(ctr_samples) if ctr_samples else "",
            ",".join(ctr_sex) if ctr_sex else "",
            ",".join(sample_samples) if sample_samples else "",
            ",".join(sample_sex) if sample_sex else "",
        )

    df[["MERGE_SEX", "MERGE_CTR", "MERGE_CTR_SEX", "MERGE_CHILDREN", "MERGE_CHILDREN_SEX"]] = \
        df["MERGE_SAMPLES"].apply(lambda x: pd.Series(extract_samples_info(x)))

    ################ classify overall SV occurance.
    df['Class'] = np.where(df['MERGE_SAMPLES'].apply(has_control), 'Ctr', 'SPonly')
    df['Class'] = np.where(df['MERGE_SAMPLES'].apply(has_child) & (df['Class'] == 'SPonly'),
                           'SPonly_children', df['Class'])

    ################ classify homo SVs.
    def extract_1_1_info(row):
        samples = row["MERGE_SAMPLES"].split(",")
        genotypes = row["MERGE_GT"].split(",")
        selected_samples = [samples[i] for i in range(len(samples)) if genotypes[i] == "1|1"]
        return ",".join(selected_samples)

    df["MERGE_1|1_All"] = df.apply(extract_1_1_info, axis=1)

    cond_Ctr = (df["MERGE_1|1_All"] != "") & df["MERGE_1|1_All"].apply(has_control)
    cond_SPonly = (df["MERGE_1|1_All"] != "") & df["MERGE_1|1_All"].apply(has_family)
    cond_SPonly_children = cond_SPonly & df["MERGE_1|1_All"].apply(has_child)
    df["HOMO_Class"] = np.select([cond_Ctr, cond_SPonly_children, cond_SPonly],
                                 ["Ctr", "SPonly_children", "SPonly"], default="")

    def extract_11_samples_info(samples, sex_dict):
        if pd.isna(samples) or not isinstance(samples, str):
            return "", "", "", "", ""
        samples = samples.strip()
        if not samples:
            return "", "", "", "", ""
        ctr_samples, child_samples, parents = [], [], []
        for s in samples.split(","):
            s = s.strip()
            g = group_of(s)
            if g == "control":
                ctr_samples.append(s)
            elif g in CHILD_GROUPS:
                child_samples.append(s)
            elif g in PARENT_GROUPS:
                parents.append(s)
        ctr_sex   = [sex_dict.get(s, "") for s in ctr_samples]
        child_sex = [sex_dict.get(s, "") for s in child_samples]
        return (
            ",".join(ctr_samples) if ctr_samples else "",
            ",".join(ctr_sex) if ctr_sex else "",
            ",".join(child_samples) if child_samples else "",
            ",".join(child_sex) if child_sex else "",
            ",".join(parents) if parents else "",
        )

    df[["MERGE_1|1_CTR", "MERGE_1|1_CTR_SEX", "MERGE_1|1_CHILDREN", "MERGE_1|1_CHILDREN_SEX", "MERGE_1|1_PARENTS"]] = \
        df["MERGE_1|1_All"].apply(lambda x: pd.Series(extract_11_samples_info(x, sex_dict)))
    oriSV = len(df)
    mis_sex = ((df['#CHROM'] =='chrY') & (df['MERGE_SEX'].str.contains('M') == False))
    df = df[~mis_sex]
    totalSV = len(df)
    print("Writing out processed input......")
    df.to_csv(os.path.join(a.outdir, 'disco_truvari_collapsed_processed.txt.gz'),
              compression='gzip', sep='\t', index=False)

    ################# extract rare_het and homo.
    rare = df[df['Class'].str.contains("SPonly_children")]
    homo = df[(df['HOMO_Class'] == "SPonly_children") &
              (df['Class'].str.contains("SPonly_children") == False)]  # included female_homo_X when ctr is het.
    rare_hom = pd.concat([rare, homo], ignore_index=True)  # non-redundant combine.

    rare_hom['SAMPLES'] = rare_hom['MERGE_SAMPLES'].str.split(',')
    rare_hom['GTS']     = rare_hom['MERGE_GT'].str.split(',')
    rare_hom['SEXS']    = rare_hom['MERGE_SEX'].str.split(',')

    def expand(frame):
        rows = []
        for _, row in frame.iterrows():
            base_data = row.to_dict()
            for sample, gt, sex in zip(row['SAMPLES'], row['GTS'], row['SEXS']):
                new_row = base_data.copy()
                new_row['Sample'] = sample
                new_row['GT'] = gt
                new_row['Sex'] = sex
                rows.append(new_row)
        out = pd.DataFrame(rows)
        return out.drop(columns=['SAMPLES', 'GTS', 'SEXS'])

    expanded_rare_hom = expand(rare_hom)
    is_child = expanded_rare_hom['Sample'].map(lambda s: group_of(s) in CHILD_GROUPS)
    rare_hom_target = expanded_rare_hom[is_child]  # retain all children rows.
    rare_hom_target_sub = rare_hom_target[~((rare_hom_target['GT'] != '1|1') & (rare_hom_target['Class'] == 'Ctr'))]

    # correct male's X GT.
    rare_hom_target_sub = rare_hom_target_sub.copy()
    rare_hom_target_sub['GT'] = np.where(
        (rare_hom_target_sub['Sex'] == 'M') & (rare_hom_target_sub['GT'] == '1|1') & (rare_hom_target_sub['#CHROM'] == 'chrX'),
        '0|1', rare_hom_target_sub['GT'])
    rare_hom_target_sub['Sex_in_Ctr'] = np.where(
        (rare_hom_target_sub['Sex'] == 'M') & (rare_hom_target_sub['GT'] == '0|1') &
        (rare_hom_target_sub['#CHROM'] == 'chrX') & (rare_hom_target_sub['MERGE_CTR_SEX'].str.contains('M') == True),
        True, False)
    rare_hom_target_sub_corr = rare_hom_target_sub[(rare_hom_target_sub['Sex_in_Ctr'] == False)]

    ################# sex-specific retrieve by sex-matching method.
    chrX = df[(df['#CHROM'] == 'chrX') &
              (df['MERGE_SAMPLES'].apply(has_child)) &
              (df['Class'] == 'Ctr') &
              (df['HOMO_Class'] != 'SPonly_children')]
    chrX = chrX.copy()
    chrX['SAMPLES'] = chrX['MERGE_SAMPLES'].str.split(',')
    chrX['GTS']     = chrX['MERGE_GT'].str.split(',')
    chrX['SEXS']    = chrX['MERGE_SEX'].str.split(',')

    expanded_chrX = expand(chrX)
    is_child_x = expanded_chrX['Sample'].map(lambda s: group_of(s) in CHILD_GROUPS)
    expanded_chrX_target = expanded_chrX[is_child_x].copy()  # retain all children rows.
    # correct male's X GT.
    expanded_chrX_target['GT'] = np.where(
        (expanded_chrX_target['Sex'] == 'M') & (expanded_chrX_target['GT'] == '1|1'),
        '0|1', expanded_chrX_target['GT'])

    ################# filtering.
    def sp_in_ctr(row):
        if row['GT'] != '1|1':
            ctr_set = set(row['MERGE_CTR_SEX'].split(','))
        else:
            ctr_set = set(row['MERGE_1|1_CTR_SEX'].split(','))
        sp_set = set(row['Sex'])
        return sp_set.issubset(ctr_set)

    expanded_chrX_target['Sex_in_Ctr'] = expanded_chrX_target.apply(sp_in_ctr, axis=1)
    chrX_target_sub = expanded_chrX_target[(expanded_chrX_target['Sex_in_Ctr'] == False)]

    ################# rare SV pool.
    pool = pd.concat([rare_hom_target_sub_corr, chrX_target_sub], ignore_index=True)
    pool['FAMID'] = pool['Sample'].map(famid_of)
    pool['role']  = pool['Sample'].map(group_of).replace({'proband': 'pro', 'sibling': 'sib'})
    pool['GenomLoc'] = np.where(pool['#CHROM'].str.contains('chrX|chrY') == False, 'auto', pool['#CHROM'])

    conditions = [(pool['GenomLoc'] == 'auto') & (pool['GT'] == '1|1'),
                  (pool['GenomLoc'] == 'auto') & (pool['GT'] != '1|1'),
                  (pool['GenomLoc'] == 'chrX') & (pool['Sex'] == 'F') & (pool['GT'] == '1|1'),
                  (pool['GenomLoc'] == 'chrX') & (pool['Sex'] == 'F') & (pool['GT'] != '1|1'),
                  (pool['GenomLoc'] == 'chrX') & (pool['Sex'] == 'M'),
                  (pool['GenomLoc'] == 'chrY') & (pool['Sex'] == 'M')]
    values = ['SPonly_Hom', 'SPonly_Het', 'SPonly_Female_HomX', 'SPonly_Female_HetX', 'SPonly_MaleX', 'SPonly_MaleY']
    pool['Category'] = np.select(conditions, values, default='unk')
    pool = pool[pool['Category'] != 'unk']
    # if your merge is "dirty", then this filters out FP.
    misY = ((pool['GenomLoc'] =='chrY') & (pool['Sex']=='M') & (pool['GT'] == '1|1') & (pool['MERGE_CTR_SEX'].str.contains('F|M') == True)|(pool['GenomLoc'] =='chrY') & (pool['Sex']=='M') & (pool['MERGE_SEX'].str.contains('F') == True)) 
    pool = pool[~misY]

    def count_families_and_members(row, samples_col):
        merge_samples = row[samples_col]
        if not isinstance(merge_samples, str) or merge_samples == "":
            return 0
        samples = merge_samples.split(',')
        # On chrY, drop any female samples that slipped in (females have no Y).
        if row['#CHROM'] == 'chrY':
            samples = [s for s in samples if sex_dict.get(s, '') != 'F']
        families = set(famid_of(s) for s in samples)
        return len(families)

    pool["Total_N_Fam"]      = pool.apply(lambda r: count_families_and_members(r, "MERGE_SAMPLES"), axis=1)
    pool["Total_N_Fam_homo"] = pool.apply(lambda r: count_families_and_members(r, "MERGE_1|1_All"), axis=1)
    print("Writing out processed rare varaint......")
    pool.to_csv(os.path.join(a.outdir, 'disco_truvari_collapsed_processed_rareSVpool.txt'), sep='\t', index=False)

    ct = pool.groupby(['Sample', 'role', 'Sex'])['ID'].count().reset_index(name='Rare')
    ct.to_csv(os.path.join(a.outdir, 'disco_truvari_collapsed_processed_rareSVpool_ct.txt'), sep='\t', index=False)
    ct_mean = ct['Rare'].mean()
    print(f"Get {ct_mean} rare SVs per child.")

    ################ optional occurrence filter.
    if a.freq is not None:
        autohet = ((pool['GenomLoc'] == 'auto') & (pool['GT'] != '1|1') & (pool["Total_N_Fam"] <= a.freq))
        autohom = ((pool['GenomLoc'] == 'auto') & (pool['GT'] == '1|1') & (pool["Total_N_Fam_homo"] <= a.freq))
        FXhet = ((pool['GenomLoc'] == 'chrX') & (pool['GT'] != '1|1') & (pool['Sex'] == 'F') & (pool["Total_N_Fam"] <= a.freq))
        FXhom = ((pool['GenomLoc'] == 'chrX') & (pool['GT'] == '1|1') & (pool['Sex'] == 'F') & (pool["Total_N_Fam_homo"] <= a.freq))
        MXYhet = ((pool['GenomLoc'] != 'auto') & (pool['Sex'] == 'M') & (pool["Total_N_Fam"] <= a.freq))

        select = autohet | autohom | FXhet | FXhom | MXYhet
        pool['select'] = np.where(select, 'y', '')
        poolsub = pool[pool["select"] == 'y']
        poolsub.to_csv(os.path.join(a.outdir, f'disco_truvari_collapsed_processed_rareSVpool{a.freq}fam.txt'), sep='\t', index=False)
        subct = poolsub.groupby(['Sample', 'role', 'Sex'])['ID'].count().reset_index(name='Rare')
        subct.to_csv(os.path.join(a.outdir, f'disco_truvari_collapsed_processed_rareSVpool{a.freq}fam_ct.txt'), sep='\t', index=False)
        subct_mean = subct['Rare'].mean()
        print(f"Get {subct_mean} rare SVs per child when freq == {a.freq}.")
        test_strct = pool[pool['Total_N_Fam'] <= a.freq]
        testct = test_strct.groupby(['Sample', 'role', 'Sex'])['ID'].count().reset_index(name='Rare')
        testct_mean = testct['Rare'].mean()
        print(f"Get {testct_mean} singletons per child when freq == {a.freq}.")

    print(f"Finished from {totalSV} total SVs ({oriSV} originally).")

if __name__ == "__main__":
    main()