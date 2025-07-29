# #!/usr/bin/env python
# encoding: utf-8
'''
@author: Yang Sui
@contact: yangsui@uw.edu
@usage: python rareSVpool.py
@input: disco_truvari_collapsed.txt.gz
@output: disco_truvari_collapsed_processed.txt.gz, disco_truvari_collapsed_processed_rareSVpool.txt
'''

import pandas as pd
import numpy as np
import sys
import os


################ add sex info.
df = pd.read_csv('example_files/108ctr_189asd_disco_truvari_collapsed.txt.gz', sep='\t')
sex_info = pd.read_csv('ctr_sp_sex.txt', sep='\t')
sex_dict = dict(zip(sex_info["sample"], sex_info["sex"]))

def extract_samples_info(samples):
    sample_list = samples.split(",")  # Split into list
    all_sex = [sex_dict.get(s, "") for s in sample_list]
    ctr_samples = [s for s in sample_list if s.startswith(("NA", "GM", "HG"))]  # Control samples
    sample_samples = [s for s in sample_list if s.endswith(("_p1", "_s1", "_s2"))]  # Other samples
    ctr_sex = [sex_dict.get(s, "") for s in ctr_samples]
    sample_sex = [sex_dict.get(s, "") for s in sample_samples]
    return ",".join(all_sex) if all_sex else "", ",".join(ctr_samples) if ctr_samples else "", ",".join(ctr_sex) if ctr_sex else "", ",".join(sample_samples) if sample_samples else "", ",".join(sample_sex) if sample_sex else ""

df[["MERGE_SEX", "MERGE_CTR", "MERGE_CTR_SEX", "MERGE_CHILDREN", "MERGE_CHILDREN_SEX"]] = df["MERGE_SAMPLES_sim"].apply(lambda x: pd.Series(extract_samples_info(x)))

################ classify overall SVs.
df['Class'] = np.where(df['MERGE_SAMPLES_sim'].str.contains('NA|GM|HG'), 'Ctr', 'Asdonly')
df['Class'] = np.where(df['MERGE_SAMPLES_sim'].str.contains('_p|_s') & (df['Class'] == 'Asdonly'), 'Asdonly_children', df['Class'])

################ classify homo SVs.
def extract_1_1_info(row):
    samples = row["MERGE_SAMPLES_sim"].split(",")  # Split sample names
    genotypes = row["MERGE_GT"].split(",")     # Split genotypes
    callers = row["MERGE_SUPP"].split(",")   # Split caller names
    selected_samples = [samples[i] for i in range(len(samples)) if genotypes[i] == "1|1"]
    selected_callers = [callers[i] for i in range(len(samples)) if genotypes[i] == "1|1"]
    return ",".join(selected_samples), ",".join(selected_callers)

df[["MERGE_1|1_All", "MERGE_1|1_Callers"]] = df.apply(extract_1_1_info, axis=1, result_type="expand")


cond_Ctr = (df["MERGE_1|1_All"] != "") & df["MERGE_1|1_All"].str.contains("NA|GM|HG")
cond_Asdonly = (df["MERGE_1|1_All"] != "") & df["MERGE_1|1_All"].str.contains("_")
cond_Asdonly_children = cond_Asdonly & df["MERGE_1|1_All"].str.contains("_p|_s")
df["HOMO_Class"] = np.select([cond_Ctr, cond_Asdonly_children, cond_Asdonly], ["Ctr", "Asdonly_children", "Asdonly"], default="")

def extract_11_samples_info(samples, sex_dict):
    if pd.isna(samples):  # Handle NaN values
        return "", "", "", ""
    if not isinstance(samples, str):  # Ensure it's a string before calling .strip()
        return "", "", "", ""
    samples = samples.strip()
    if not samples:
        return "", "", "", ""
    ctr_samples = []
    child_samples = []
    parents = []
    for s in samples.split(","):
        s = s.strip()
        if s.startswith(("NA", "GM", "HG")):
            ctr_samples.append(s)
        elif s.endswith(("_p1", "_s1", "_s2")):
            child_samples.append(s)
        elif s.endswith(("_fa", "_mo")):
            parents.append(s)
    ctr_sex = [sex_dict.get(s, "") for s in ctr_samples]
    child_sex = [sex_dict.get(s, "") for s in child_samples]
    return (
        ",".join(ctr_samples) if ctr_samples else "",
        ",".join(ctr_sex) if ctr_sex else "",
        ",".join(child_samples) if child_samples else "",
        ",".join(child_sex) if child_sex else "",
        ",".join(parents) if parents else ""
    )

df[["MERGE_1|1_CTR", "MERGE_1|1_CTR_SEX", "MERGE_1|1_CHILDREN", "MERGE_1|1_CHILDREN_SEX", "MERGE_1|1_PARENTS"]] = df["MERGE_1|1_All"].apply(lambda x: pd.Series(extract_11_samples_info(x, sex_dict)))
simdf = df[['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'MERGE_SAMPLES_sim', 'MERGE_GT', "MERGE_SEX", 
       'Class', 'MERGE_CTR', 'MERGE_CTR_SEX', 'MERGE_CHILDREN', 'MERGE_CHILDREN_SEX',
       'HOMO_Class', 'MERGE_1|1_CTR', 'MERGE_1|1_CTR_SEX', "MERGE_1|1_CHILDREN", "MERGE_1|1_CHILDREN_SEX", "MERGE_1|1_PARENTS", 'MERGE_1|1_All']]
# simdf.to_csv('disco_truvari_collapsed_processed.txt.gz', compression='gzip', sep='\t', index=False)



################# extract rare_het and homo.
rare = simdf[simdf['Class'].str.contains("Asdonly_children")]
homo = simdf[(simdf['HOMO_Class'] == "Asdonly_children") & (simdf['Class'].str.contains("Asdonly_children") == False)] # included female_homo_X when ctr is het.
rare_hom = pd.concat([rare, homo], ignore_index=True) # it's a non-redundant combine, so "rare" might include some "homo"; and homo might include maleX_homo cases.

rare_hom['SAMPLES'] = rare_hom['MERGE_SAMPLES_sim'].str.split(',')
rare_hom['GTS'] = rare_hom['MERGE_GT'].str.split(',')
rare_hom['SEXS'] = rare_hom['MERGE_SEX'].str.split(',')
expanded_rows = []

for _, row in rare_hom.iterrows():
    base_data = row.to_dict()  # Convert entire row to dictionary
    merge_samples = row['SAMPLES']
    merge_gts = row['GTS']
    merge_sexs = row['SEXS']
    
    # Expand for each Sample-GT pair
    for sample, gt, sex in zip(merge_samples, merge_gts, merge_sexs):
        new_row = base_data.copy()  # Retain all original columns
        new_row['Sample'] = sample
        new_row['GT'] = gt
        new_row['Sex'] = sex
        expanded_rows.append(new_row)


expanded_rare_hom = pd.DataFrame(expanded_rows)
expanded_rare_hom.drop(columns=['SAMPLES', 'GTS', 'SEXS'], inplace=True)
rare_hom_target = expanded_rare_hom[expanded_rare_hom['Sample'].str.contains('_p|_s')] # Retain all children rows.
rare_hom_target_sub = rare_hom_target[~((rare_hom_target['GT'] != '1|1') & (rare_hom_target['Class'] == 'Ctr'))] 
# correct male's X GT.
rare_hom_target_sub['GT'] = np.where((rare_hom_target_sub['Sex'] == 'M') & (rare_hom_target_sub['GT'] == '1|1') & (rare_hom_target_sub['#CHROM'] == 'chrX'), '0|1', rare_hom_target_sub['GT'])
rare_hom_target_sub['Sex_in_Ctr'] = np.where((rare_hom_target_sub['Sex'] == 'M') & (rare_hom_target_sub['GT'] == '0|1') & (rare_hom_target_sub['#CHROM'] == 'chrX') & (rare_hom_target_sub['MERGE_CTR_SEX'].str.contains('M') == True), True, False)
rare_hom_target_sub_corr = rare_hom_target_sub[(rare_hom_target_sub['Sex_in_Ctr'] == False)] 

################# sex-specific retrieve by sex-matching method.
chrX = simdf[(simdf['#CHROM'] == 'chrX') & (simdf['MERGE_SAMPLES_sim'].str.contains('_p|_s')) & (simdf['Class'] == 'Ctr') & (simdf['HOMO_Class'] != 'Asdonly_children')]
chrX['SAMPLES'] = chrX['MERGE_SAMPLES_sim'].str.split(',')
chrX['GTS'] = chrX['MERGE_GT'].str.split(',')
chrX['SEXS'] = chrX['MERGE_SEX'].str.split(',')
expanded_chrX = []

for _, row in chrX.iterrows():
    base_data = row.to_dict()  # Convert entire row to dictionary
    merge_samples = row['SAMPLES']
    merge_gts = row['GTS']
    merge_sexs = row['SEXS']
    
    # Expand for each Sample-GT pair
    for sample, gt, sex in zip(merge_samples, merge_gts, merge_sexs):
        new_row = base_data.copy()  # Retain all original columns
        new_row['Sample'] = sample
        new_row['GT'] = gt
        new_row['Sex'] = sex
        expanded_chrX.append(new_row)


expanded_chrX = pd.DataFrame(expanded_chrX)
expanded_chrX.drop(columns=['SAMPLES', 'GTS', 'SEXS'], inplace=True)
expanded_chrX_target = expanded_chrX[expanded_chrX['Sample'].str.contains('_p|_s')] # Retain all children rows.
# correct male's X GT.
expanded_chrX_target['GT'] = np.where((expanded_chrX_target['Sex'] == 'M') & (expanded_chrX_target['GT'] == '1|1'), '0|1', expanded_chrX_target['GT'])

################# filtering.
def sp_in_ctr(row):
    if row['GT'] != '1|1':
        ctr_set = set(row['MERGE_CTR_SEX'].split(','))
    else:  # GT == '1|1'
        ctr_set = set(row['MERGE_1|1_CTR_SEX'].split(','))
    sp_set = set(row['Sex'])  # assumes 'Sex' is a string like 'F' or 'M'
    return sp_set.issubset(ctr_set)

expanded_chrX_target['Sex_in_Ctr'] = expanded_chrX_target.apply(sp_in_ctr, axis=1)
chrX_target_sub = expanded_chrX_target[(expanded_chrX_target['Sex_in_Ctr'] == False)] 


################# rare SV pool。
pool= pd.concat([rare_hom_target_sub_corr, chrX_target_sub], ignore_index=True) 
pool['FAMID'] = pool['Sample'].str.split('_').str[0]
pool['role'] = pool['Sample'].str.split('_').str[1]
pool['role'] = pool['role'].str.replace('p1', 'pro').str.replace('s1', 'sib').str.replace('s2', 'sib')
pool['GenomLoc'] = np.where(pool['#CHROM'].str.contains('chrX|chrY') == False, 'auto', pool['#CHROM'])

################# add parental GTs, correct children's and fatherX's GTs。
def extract_family_id(sample):
    return sample.split("_")[0]  # Extracts 'aa' or 'bb'

def extract_parent_GT(row):
    samples = row["MERGE_SAMPLES_sim"].split(",")
    GTs = row["MERGE_GT"].split(",")
    sample_gt_map = dict(zip(samples, GTs))
    family_id = extract_family_id(row["Sample"])
    faGT = sample_gt_map.get([s for s in samples if s.startswith(family_id) and "_fa" in s][0], "") if any(s.startswith(family_id) and "_fa" in s for s in samples) else ""
    moGT = sample_gt_map.get([s for s in samples if s.startswith(family_id) and "_mo" in s][0], "") if any(s.startswith(family_id) and "_mo" in s for s in samples) else ""
    return pd.Series([faGT, moGT])


pool[["faGT_HC", "moGT_HC"]] = pool.apply(extract_parent_GT, axis=1).replace("", np.nan)
pool['GT'] = np.where((pool['#CHROM'] == 'chrX') & (pool['Sex'] == 'M') & (pool['GT'] != '0|1'), '0|1', pool['GT']) # correct maleX GTs.
pool['GT'] = np.where((pool['#CHROM'] == 'chrY') & (pool['Sex'] == 'M') & (pool['GT'] != '1|0'), '1|0', pool['GT']) # correct maleY GTs.
pool['faGT_HC'] = np.where((pool['#CHROM'] == 'chrX') & (pool['faGT_HC'] == '1|1'), '0|1', pool['faGT_HC'])
pool['faGT_HC'] = np.where((pool['#CHROM'] == 'chrY') & (pool['faGT_HC'].str.contains('1') == True), '1|0', pool['faGT_HC'])


################# add category.
conditions = [(pool['GenomLoc'] == 'auto') & (pool['GT'] == '1|1'),
              (pool['GenomLoc'] == 'auto') & (pool['GT'] != '1|1'),
              (pool['GenomLoc'] == 'chrX') & (pool['Sex'] == 'F') & (pool['GT'] == '1|1'),
              (pool['GenomLoc'] == 'chrX') & (pool['Sex'] == 'F') & (pool['GT'] != '1|1'),
              (pool['GenomLoc'] == 'chrX') & (pool['Sex'] == 'M'),
              (pool['GenomLoc'] == 'chrY') & (pool['Sex'] == 'M')]

values = ['Asdonly_Hom', 'Asdonly_Het', 'Asdonly_Female_HomX', 'Asdonly_Female_HetX','Asdonly_MaleX', 'Asdonly_MaleY']
pool['Category'] = np.select(conditions, values, default='unk')

################# add family type.
famid = pd.read_csv('famid.txt', sep='\t')
pool = pd.merge(pool, famid, on='FAMID', how='left')
print(len(pool))
def extract_famps(row):
    samples = row['MERGE_SAMPLES_sim'].split(',')
    gts = row['MERGE_GT'].split(',')
    famid = row['FAMID']
    pro_gt = '0|0'
    sib_gt = '0|0'   
    for sample, gt in zip(samples, gts):
        if sample.startswith(famid):
            if '_p' in sample:
                pro_gt = gt
            elif '_s' in sample:
                sib_gt = gt    
    return pd.Series({'pro_GT': pro_gt, 'sib_GT': sib_gt})

pool[['pro_GT', 'sib_GT']] = pool.apply(extract_famps, axis=1)
pool['faGT_HC'] = pool['faGT_HC'].fillna('0|0')
pool['moGT_HC'] = pool['moGT_HC'].fillna('0|0')

def check_mendelian(row):
    fa = set(row['faGT_HC'].split('|'))
    mo = set(row['moGT_HC'].split('|'))
    possible_alleles = {f + '|' + m for f in fa for m in mo}
    def identity(gt):
        return gt
    gt = identity(row['GT'])
    if gt in possible_alleles:
        return 'HC'
    else:
        return 'LC_MENDEL'

pool['Conf'] = pool.apply(check_mendelian, axis=1)


################ count frequency.
def count_families_and_members(merge_samples):
    if merge_samples == "" or not isinstance(merge_samples, str):
        return 0  # Handle empty string or non-string cases gracefully
    families = set(sample.split('_')[0] for sample in merge_samples.split(','))
    return len(families)

pool["Total_N_Fam"] = pool["MERGE_SAMPLES_sim"].apply(count_families_and_members)
pool["Total_N_Fam_homo"] = pool["MERGE_1|1_All"].apply(count_families_and_members)

HC_Singleton = (pool['Conf'] == 'HC') & (pool['Total_N_Fam'] == 1)

private = [
(pool['Category'] == 'Asdonly_MaleX') & (pool['moGT_HC'].str.contains('0')) & (pool['faGT_HC'] == '0|0') & (pool['Class'] != 'Ctr') & HC_Singleton,
(pool['Category'].str.contains('Het') == True) & HC_Singleton & (pool['moGT_HC'].str.contains('0') == True) & (pool['faGT_HC'] == '0|0') & (pool['GT'] == '0|1'),
(pool['Category'].str.contains('Het') == True) & HC_Singleton & (pool['faGT_HC'].str.contains('0') == True) & (pool['moGT_HC'] == '0|0') & (pool['GT'] == '1|0')]
pi_value = ['PISV', 'PISV', 'PISV']
pool['Private'] = np.select(private, pi_value, default='')
pool['Private'] = np.where(HC_Singleton & (pool['GenomLoc'] == 'chrY'), 'PISV', pool['Private']) # add private chrY.
pool['Private'] = np.where((pool["Total_N_Fam_homo"] ==1) & (pool['GT'] == '1|1') & (pool['Conf'] == 'HC') & (pool['MERGE_1|1_PARENTS']  == ''), 'PISV_hom', pool['Private']) # add privtae hom.

pool.to_csv('example_files/disco_truvari_collapsed_processed_rareSVpool.txt', sep='\t', index=False)










