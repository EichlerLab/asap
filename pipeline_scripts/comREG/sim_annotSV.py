# #!/usr/bin/env python

# encoding: utf-8

'''
@author: Yang Sui
@contact: yangsui@uw.edu
@ Prerequisite: 
module load annotsv/3.4
AnnotSV -SvinputFile annotSV.input.bed -svtBEDcol 4 -candidateGenesFile /net/eichler/vol28/projects/autism_genome_assembly/nobackups/yangsui/annotSV/NDDgene.txt -promoterSize 500 -REreport 1 -outputFile truvariSVsanno
'''
import pandas as pd
import numpy as np

# read file by parsing lines.
finp = open('example_files/truvariSVsanno.tsv')
header = finp.readline().strip().split('\t')
data = [ line.replace('\n','').split('\t') for line in finp ]
df = pd.DataFrame(data)
df.columns=header

print(df.info())

sub = df[["user#1","Annotation_mode","Gene_name","Gene_count","Location","SV_chrom","SV_start","SV_end","CytoBand","Closest_left","Closest_right","RE_gene","P_gain_phen","P_loss_phen","P_ins_phen", "po_P_gain_phen", "po_P_loss_phen","GC_content_left","GC_content_right","Repeat_type_left","Repeat_type_right","ACMG","HI","TS","DDD_HI_percent","GenCC_disease","GenCC_moi","GenCC_classification","GenCC_pmid","NCBI_gene_ID","OMIM_ID","OMIM_phenotype","OMIM_inheritance","OMIM_morbid","OMIM_morbid_candidate","LOEUF_bin","GnomAD_pLI","ExAC_pLI","AnnotSV_ranking_score","ACMG_class"]]

full = sub[sub['Annotation_mode'] == 'full']
uniqfull = full.drop_duplicates(keep='first')

split = sub[sub['Annotation_mode'] == 'split']
uniqsplit = split.drop_duplicates(keep='first')

uniqsplit['HI-TS'] = uniqsplit['HI'].astype(str) + '_' + uniqsplit['TS'].astype(str)
uniqsplit['GenCC'] = uniqsplit['GenCC_disease'].astype(str) + '|' + uniqsplit['GenCC_moi'].astype(str) + '|' + uniqsplit['GenCC_classification'].astype(str) + '|' + uniqsplit['GenCC_pmid'].astype(str)
uniqsplit['OMIM'] = uniqsplit['OMIM_phenotype'].astype(str) + '|' + uniqsplit['OMIM_inheritance'].astype(str) + '|' + uniqsplit['OMIM_morbid'].astype(str) + '|' + uniqsplit['OMIM_morbid_candidate'].astype(str)

uniqsplit['left'] = uniqsplit['Location'].str.split('-').str[0]
uniqsplit['right'] = uniqsplit['Location'].str.split('-').str[1]
uniqsplit['LocationSV'] = np.where((uniqsplit['left'] == uniqsplit['right']) & (uniqsplit['left'].str.contains('intron')), 'Intronic', 'Exonic')


sim_uniqsplit = uniqsplit.groupby("user#1", as_index=False).agg({
    "Gene_name": ";".join,
    "LocationSV": ";".join,
    "HI-TS": "//".join,
    "GenCC": "//".join,
    "OMIM": "//".join
})

merge = pd.merge(uniqfull, sim_uniqsplit, on='user#1', how='left')
print(len(merge),len(uniqfull))

merge['geneMatch'] = np.where(merge['Gene_name_x'] == merge['Gene_name_y'], 1, 0)
merge['Gene_count']= merge['Gene_count'].astype(int)
merge['LocationSV'] = merge['LocationSV'].astype(str)
conditions = [(merge['Gene_count'] == 0), (merge['LocationSV'].str.contains('Exonic', na=False)), (merge['LocationSV'].str.contains('Intronic', na=False))]
values = ['Intergenic', 'Exonic', 'Intronic']
merge['LocationSV_sim'] = np.select(conditions, values)


final_anno = merge[['user#1', 'LocationSV_sim', 'geneMatch', 'Gene_count','Gene_name_x', 'Gene_name_y', 'LocationSV', 'AnnotSV_ranking_score', 
       'ACMG_class', 'HI-TS', 'GenCC', 'OMIM', 'CytoBand', 'GC_content_left', 'GC_content_right', 'Closest_left', 'Closest_right', 'RE_gene', 'P_gain_phen',
       'P_loss_phen', 'P_ins_phen', 'po_P_gain_phen', 'po_P_loss_phen','Repeat_type_left', 'Repeat_type_right', 'ACMG', 'HI', 'TS', 'DDD_HI_percent', 'LOEUF_bin', 'GnomAD_pLI','ExAC_pLI']]


final_anno = final_anno.rename({'user#1': 'ID'}, axis=1)

final_anno.to_csv('example_files/truvariSVsanno_sim.txt', sep='\t', index=False)











