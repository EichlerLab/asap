import pandas as pd
import numpy as np

sv_file = 'annotation_reg_input.bed'
manifest_df = pd.read_csv('manifest.txt', sep='\t', header=0, index_col="peak")

# pli = '/net/eichler/vol28/projects/autism_genome_assembly/nobackups/yangsui/annotation/data/Gnomad_constraint/largest_pLI_coor.txt'
# loeuf = '/net/eichler/vol28/projects/autism_genome_assembly/nobackups/yangsui/annotation/data/Gnomad_constraint/largest_loeuf_coor.txt'
# ndd = '/net/eichler/vol28/projects/autism_genome_assembly/nobackups/yangsui/annotation/data/ASD/anno2.xlsx'
# hpo = '/net/eichler/vol28/projects/autism_genome_assembly/nobackups/yangsui/annotation/data/HPO/HPOsim.txt'
# pheno = '/net/eichler/vol28/projects/autism_genome_assembly/nobackups/yangsui/annotation/data/HPO/Pheno_ENSEMBL.txt'

wildcard_constraints:
    peak='|'.join(manifest_df.index)

def get_peak(wildcards):
    return manifest_df.at[wildcards.peak, "path"]

def parse_gencode(gencode):
    if pd.isna(gencode) or gencode == '':
        return '', '', '', ''
    all_genes = []
    cds, utr5, utr3 = [], [], []    
    for item in gencode.split(';'):
        if '_' not in item:
            continue
        region, gene = item.split('_', 1)        
        # all genes (unique, preserve order)
        if gene not in all_genes:
            all_genes.append(gene)        
        # region-specific
        if region == 'CDS' and gene not in cds:
            cds.append(gene)
        elif region == '5UTR' and gene not in utr5:
            utr5.append(gene)
        elif region == '3UTR' and gene not in utr3:
            utr3.append(gene)    
    return (
        ';'.join(all_genes),
        ';'.join(cds),
        ';'.join(utr5),
        ';'.join(utr3),
    )

def get_gene_values(gene, values_df, column_name):
    if pd.isna(gene):  # Check for NaN values
        return pd.NA
    if ';' in gene:
        result = ';'.join(
            f"{g}:{values_df.loc[values_df['gene'] == g, column_name].values[0]}"
            for g in gene.split(';') if g in values_df['gene'].values
        )
        return result if result else pd.NA  # Return NA if result is empty
    else:
        return f"{gene}:{values_df.loc[values_df['gene'] == gene, column_name].values[0]}" \
            if gene in values_df['gene'].values else pd.NA

def combine_or_no(row):
    values = [str(x) for x in row if pd.notna(x) and str(x).strip() != '']
    return ';'.join(values) if values else '.'


rule all:
    input:
        expand(['allreg_sv.txt'], peak=manifest_df.index)

#
# variant_anno_intersect
#
# Intersect variants with PEAKS.
rule variant_anno_intersect:
    input:
        bed=sv_file,
        anno=get_peak,
    output:
        tsv='temp/{peak}_sv.txt',
    resources:
        mem=10,
        hrs=24,
    envmodules:
        "miniconda/4.12.0"
    threads: 1
    #benchmark: "benchmark/oreganno.log"
    shell:
        """
        bedtools intersect -a {input.bed} -b {input.anno} -loj| awk -v OFS='\t' '{{ if($5 != ".") print $0;}}' |cut -f4,8 | sort | uniq | awk '{{a[$1]=a[$1] ? a[$1] ";" $2 : $2}} END {{for (i in a) print i"\t"a[i]}}' > {output.tsv}
        """


rule TRSV:
    input:
        bed=sv_file,
        TR=manifest_df.at["platinumTRs", "path"],
    output:
        tsv='temp/TRsv.txt',
    resources:
        mem=10,
        hrs=24,
    envmodules:
        "miniconda/4.12.0"
    threads: 1
    #benchmark: "benchmark/oreganno.log"
    shell:
        """
        bedtools intersect -a {input.bed} -b {input.TR} -f 0.5 -wa -wb| cut -f4,8 | sort | uniq | awk '{{a[$1]=a[$1] ? a[$1] ";" $2 : $2}} END {{for (i in a) print i"\t"a[i]}}' > {output.tsv}
        """


rule variant_anno_max_constraint:
    input:
        bed=sv_file,
        nc_constraint='data/GRCh38/Non-coding_constraint_Supplementary_Data_2.bed'
    output:
        nc_constraint_tsv='temp/nc_constraint_sv.txt'
    resources:
        mem=10,
        hrs=24,
    envmodules:
        "miniconda/4.12.0"
    threads: 1
    #benchmark: "benchmark/oreganno.log"
    shell:
        """
        bedtools intersect -a {input.bed} -b {input.nc_constraint} -wo|sort -k4,4 -k8,8r|awk '!seen[$4]++' |cut -f4,8 > {output.nc_constraint_tsv}
        """

rule variant_anno_max_segdup:
    input:
        bed=sv_file,
        segdup='data/GRCh38/SegDup.bed'
    output:
        segdup_tsv='temp/segdup_sv.txt'
    resources:
        mem=10,
        hrs=24,
    envmodules:
        "miniconda/4.12.0"
    threads: 1
    #benchmark: "benchmark/oreganno.log"
    shell:
        """
        bedtools intersect -a {input.bed} -b {input.segdup} -wo|sort -k4,4 -k8,8r|awk '!seen[$4]++' |cut -f4,8 > {output.segdup_tsv}
        """

#
# All regulatory
#

# variant_anno_reg_all_reg
#
# Merge all regulatory elements into one table.
rule variant_anno_reg_all_reg:
    input:
        bed=sv_file,
        ccre='temp/cCREs_sv.txt',
        genehancer='temp/geneHancer_sv.txt',
        oreganno='temp/oreganno_sv.txt',
        atac='temp/CorticalMap_sv.txt',
        H3K27Ac='temp/H3K27Ac_sv.txt',
        H3K4Me1='temp/H3K4Me1_sv.txt',
        H3K4Me3='temp/H3K4Me3_sv.txt',
        TF='temp/TFcluster_sv.txt',
        K_H3K27ac='temp/K_H3K27ac_sv.txt',
        K_H3K27me3='temp/K_H3K27me3_sv.txt',
        K_H3K4me3='temp/K_H3K4me3_sv.txt',
        K_CTCF='temp/K_CTCF_sv.txt',
        TRSV='temp/TRsv.txt',
        platinumTRs='temp/platinumTRs_sv.txt',
        RepeatMasker='temp/RepeatMasker_sv.txt',
        SimpleRepeat='temp/SimpleRepeat_sv.txt',
        tRNA='temp/tRNA_sv.txt',
        snRNA='temp/snRNA_sv.txt',
        lincRNA='temp/lincRNA_sv.txt',
        sno_miRNA='temp/sno_miRNA_sv.txt',
        segdup=rules.variant_anno_max_segdup.output.segdup_tsv,
        Gnocchi=rules.variant_anno_max_constraint.output.nc_constraint_tsv,
        gene = 'temp/gene_sv.txt',
        CpGI = 'temp/CpGI_sv.txt',
        morbiditymapDEL = 'temp/mDEL_sv.txt',
        morbiditymapDUP = 'temp/mDUP_sv.txt',
        # pli = pli,
        # loeuf = loeuf,
        # ndd = ndd,
        # hpo = hpo,
        # pheno = pheno
    output:
        tsv='allreg_sv.txt',
    resources:
        mem=20,
        hrs=24,
    threads: 1,
    #benchmark: "benchmark/all.log"
    run:
        df = pd.read_csv(input.bed, sep='\t')
        annotation_columns = ['GENCODE', 'Gnocchi', 'ENCODE_cCRE', 'ENCODE_H3K27Ac', 'ENCODE_H3K4Me1', 'ENCODE_H3K4Me3', 'ORegAnno', 'GeneHancer', 'ENCODE_TFcluster', 'brain_ATAC_CorticalMap', 'brain_CUT&Tag_H3K27ac', 'brain_CUT&Tag_H3K27me3', 'brain_CUT&Tag_H3K4me3', 'brain_CUT&Tag_CTCF', 'SegDup', 'RepeatMasker', 'SimpleRepeat', 'platinumTRs', 'TRSV','tRNA', 'snRNA', 'lincRNA', 'sno_miRNA','CpGI','morbiditymapDEL','morbiditymapDUP']
        for col in annotation_columns:
            df[col] = None  # Initialize all columns with None

        df.set_index('ID', inplace=True)
        annotation_files = {
            'GENCODE':input.gene,
            'ENCODE_cCRE':input.ccre,
            'GeneHancer':input.genehancer,
            'ORegAnno':input.oreganno,
            'ENCODE_H3K27Ac':input.H3K27Ac,
            'ENCODE_H3K4Me1':input.H3K4Me1,
            'ENCODE_H3K4Me3':input.H3K4Me3,
            'ENCODE_TFcluster':input.TF,
            'brain_ATAC_CorticalMap':input.atac,
            'brain_CUT&Tag_H3K27ac':input.K_H3K27ac,
            'brain_CUT&Tag_H3K27me3':input.K_H3K27me3,
            'brain_CUT&Tag_H3K4me3':input.K_H3K4me3,
            'brain_CUT&Tag_CTCF':input.K_CTCF,
            'platinumTRs':input.platinumTRs,
            'TRSV':input.TRSV,
            'RepeatMasker':input.RepeatMasker,
            'SimpleRepeat':input.SimpleRepeat,
            'tRNA':input.tRNA,
            'snRNA':input.snRNA,
            'lincRNA':input.lincRNA,
            'sno_miRNA':input.sno_miRNA,
            'SegDup':input.segdup,
            'Gnocchi':input.Gnocchi,
            'CpGI':input.CpGI,
            'morbiditymapDEL':input.morbiditymapDEL,
            'morbiditymapDUP':input.morbiditymapDUP
            }

        for column, file_path in annotation_files.items():
            annotation_df = pd.read_csv(file_path, sep='\t', header=None, names=['ID', 'value'])
            annotation_df.set_index('ID', inplace=True)
            df.loc[df.index.isin(annotation_df.index), column] = annotation_df.loc[df.index.intersection(annotation_df.index), 'value']

        df['REG'] = df[['ENCODE_cCRE','ENCODE_H3K27Ac', 'ENCODE_H3K4Me1', 'ENCODE_H3K4Me3','ORegAnno', 'GeneHancer', 'brain_ATAC_CorticalMap','brain_CUT&Tag_H3K27ac', 'brain_CUT&Tag_H3K27me3', 'brain_CUT&Tag_H3K4me3','brain_CUT&Tag_CTCF']].notna().any(axis=1).replace({True: 'YES', False: 'NO'})
        df['REGwTF'] = df[['ENCODE_cCRE','ENCODE_H3K27Ac', 'ENCODE_H3K4Me1', 'ENCODE_H3K4Me3','ORegAnno', 'GeneHancer', 'brain_ATAC_CorticalMap','brain_CUT&Tag_H3K27ac', 'brain_CUT&Tag_H3K27me3', 'brain_CUT&Tag_H3K4me3','brain_CUT&Tag_CTCF','ENCODE_TFcluster']].notna().any(axis=1).replace({True: 'YES', False: 'NO'})
        df['otherREG'] = df[['ENCODE_cCRE','ENCODE_H3K27Ac', 'ENCODE_H3K4Me1', 'ENCODE_H3K4Me3','ORegAnno', 'GeneHancer']].notna().any(axis=1).replace({True: 'YES', False: 'NO'})
        df['brainREG'] = df[['brain_ATAC_CorticalMap','brain_CUT&Tag_H3K27ac', 'brain_CUT&Tag_H3K27me3', 'brain_CUT&Tag_H3K4me3','brain_CUT&Tag_CTCF']].notna().any(axis=1).replace({True: 'YES', False: 'NO'})
        df[['gene', 'CDS_genes', '5UTR_genes', '3UTR_genes']] = (df['GENCODE'].apply(parse_gencode).apply(pd.Series))
        anno_reg_annot = df
        # pli = pd.read_csv(input.pli, sep='\t')
        # pli['lof.pLI_gnomad.v4.1'] = pli['lof.pLI_gnomad.v4.1'].round(3)
        # loeuf = pd.read_csv(input.loeuf, sep='\t')
        # loeuf['lof.oe_ci.upper_gnomad.v4.1'] = loeuf['lof.oe_ci.upper_gnomad.v4.1'].round(3)
        # hpo = pd.read_csv(input.hpo, sep='\t')
        # pheno = pd.read_csv(input.pheno, sep='\t')
        # sfari = pd.read_excel(input.ndd, sheet_name='SFARI07082025')
        # fu = pd.read_excel(input.ndd, sheet_name='NDD_Fu')
        # sat = pd.read_excel(input.ndd, sheet_name='ASD102')
        # tw = pd.read_excel(input.ndd, sheet_name='NDD_tw')
        # ys = pd.read_excel(input.ndd, sheet_name='ys')
        # mrg = pd.read_excel(input.ndd, sheet_name='MRG')
        # cmrg = pd.read_excel(input.ndd, sheet_name='373CMRG')
        # gene_mappings = {'pLIv4.1': (pli, 'lof.pLI_gnomad.v4.1'),'LOEUFv4.1': (loeuf, 'lof.oe_ci.upper_gnomad.v4.1'),'SFARI_070825': (sfari, 'SFARI_070825'),'Fu_NDD664': (fu, 'Fu_NDD664'),'Wang_NDD615': (tw, 'tw_NDD615'),'Satterstrom_ASD102': (sat, 'Satt_ASD102'),'Sui_NDD810': (ys, 'Sui_810'),'MRG': (mrg, 'MRG'), 'CMRG': (cmrg, 'CMRG'), 'HPO_name': (hpo, 'hpo_name'), 'Disease_ID': (hpo, 'disease_id'), 'Phenotype': (pheno, 'Phenotype')}
        # for new_col, (df, col_name) in gene_mappings.items():
        #     anno_reg_annot[new_col] = anno_reg_annot['gene'].apply(get_gene_values, args=(df, col_name))

        # anno_reg_annot['NDDsui'] = anno_reg_annot[['SFARI_070825', 'Fu_NDD664', 'Satterstrom_ASD102', 'Wang_NDD615', 'Sui_NDD810']].apply(combine_or_no, axis=1)
        # anno_reg_annot['NDD'] = anno_reg_annot[['SFARI_070825', 'Fu_NDD664', 'Satterstrom_ASD102', 'Wang_NDD615']].apply(combine_or_no, axis=1)
        anno_reg_annot['Location_GENCODE'] = '.'
        anno_reg_annot.loc[anno_reg_annot['3UTR_genes'] != '', 'Location_GENCODE'] = '3UTR'
        anno_reg_annot.loc[anno_reg_annot['5UTR_genes'] != '', 'Location_GENCODE'] = '5UTR'
        anno_reg_annot.loc[anno_reg_annot['CDS_genes']  != '', 'Location_GENCODE'] = 'CDS'
        anno_reg_annot['Location_GENCODE'] = np.where((anno_reg_annot['Location_GENCODE'] == '.') & (anno_reg_annot['REG'] == 'YES'), 'REG', anno_reg_annot['Location_GENCODE'])
        anno_reg_annot.to_csv(output.tsv, sep='\t', index=True)


