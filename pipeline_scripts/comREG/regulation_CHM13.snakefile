import pandas as pd
import numpy as np

sv_file = 'annotation_reg_input.bed'
manifest_df = pd.read_csv('manifestCHM13.txt', sep='\t', header=0, index_col="peak")

# pli = '/net/eichler/vol28/projects/autism_genome_assembly/nobackups/yangsui/annotation/data/Gnomad_constraint/largest_pLI_coor.txt'
# loeuf = '/net/eichler/vol28/projects/autism_genome_assembly/nobackups/yangsui/annotation/data/Gnomad_constraint/largest_loeuf_coor.txt'
# ndd = '/net/eichler/vol28/projects/autism_genome_assembly/nobackups/yangsui/annotation/data/ASD/anno2.xlsx'
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
        expand(['allreg_sv_T2T.txt'], peak=manifest_df.index)

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

# All regulatory
#

# variant_anno_reg_all_reg
#
# Merge all regulatory elements into one table.
rule variant_anno_reg_all_reg:
    input:
        bed=sv_file,
        ccre='temp/cCREs_sv.txt',
        platinumTRs='temp/TRsv.txt',
        RepeatMasker='temp/RepeatMasker_sv.txt',
        SimpleRepeat='temp/SimpleRepeat_sv.txt',
        segdup='temp/SD_sv.txt',
        refseq='temp/gene_sv.txt',
        CpGI = 'temp/CpGI_sv.txt',
        morbiditymapDEL = 'temp/mDEL_sv.txt',
        morbiditymapDUP = 'temp/mDUP_sv.txt',
        # pli = pli,
        # loeuf = loeuf,
        # ndd = ndd
    output:
        tsv='allreg_sv_T2T.txt',
    resources:
        mem=20,
        hrs=24,
    threads: 1,
    #benchmark: "benchmark/all.log"
    run:
        df = pd.read_csv(input.bed, sep='\t')
        annotation_columns = ['ENCODE_cCRE', 'SegDup', 'RepeatMasker', 'SimpleRepeat', 'platinumTRs','refseq', 'CpGI','morbiditymapDEL','morbiditymapDUP']
        for col in annotation_columns:
            df[col] = None  # Initialize all columns with None

        df.set_index('ID', inplace=True)
        annotation_files = {
            'ENCODE_cCRE':input.ccre,
            'SegDup':input.segdup,
            'RepeatMasker':input.RepeatMasker,
            'SimpleRepeat':input.SimpleRepeat,
            'platinumTRs':input.platinumTRs,
            'refseq':input.refseq,
            'CpGI':input.CpGI,
            'morbiditymapDEL':input.morbiditymapDEL,
            'morbiditymapDUP':input.morbiditymapDUP
            }

        for column, file_path in annotation_files.items():
            annotation_df = pd.read_csv(file_path, sep='\t', header=None, names=['ID', 'value'])
            annotation_df.set_index('ID', inplace=True)
            df.loc[df.index.isin(annotation_df.index), column] = annotation_df.loc[df.index.intersection(annotation_df.index), 'value']

        df[['gene', 'CDS_genes', '5UTR_genes', '3UTR_genes']] = (df['refseq'].apply(parse_gencode).apply(pd.Series))
        anno_reg_annot = df
        # pli = pd.read_csv(input.pli, sep='\t')
        # pli['lof.pLI_gnomad.v4.1'] = pli['lof.pLI_gnomad.v4.1'].round(3)
        # loeuf = pd.read_csv(input.loeuf, sep='\t')
        # loeuf['lof.oe_ci.upper_gnomad.v4.1'] = loeuf['lof.oe_ci.upper_gnomad.v4.1'].round(3)
        # sfari = pd.read_excel(input.ndd, sheet_name='SFARI07082025')
        # fu = pd.read_excel(input.ndd, sheet_name='NDD_Fu')
        # sat = pd.read_excel(input.ndd, sheet_name='ASD102')
        # tw = pd.read_excel(input.ndd, sheet_name='NDD_tw')
        # ys = pd.read_excel(input.ndd, sheet_name='ys')
        # gene_mappings = {'pLIv4.1': (pli, 'lof.pLI_gnomad.v4.1'),'LOEUFv4.1': (loeuf, 'lof.oe_ci.upper_gnomad.v4.1'),'SFARI_070825': (sfari, 'SFARI_070825'),'Fu_NDD664': (fu, 'Fu_NDD664'),'Wang_NDD615': (tw, 'tw_NDD615'),'Satterstrom_ASD102': (sat, 'Satt_ASD102'),'Sui_NDD810': (ys, 'Sui_810')}
        # for new_col, (df, col_name) in gene_mappings.items():
        #     anno_reg_annot[new_col] = anno_reg_annot['gene'].apply(get_gene_values, args=(df, col_name))

        # anno_reg_annot['NDDsui'] = anno_reg_annot[['SFARI_070825', 'Fu_NDD664', 'Satterstrom_ASD102', 'Wang_NDD615', 'Sui_NDD810']].apply(combine_or_no, axis=1)
        # anno_reg_annot['NDD'] = anno_reg_annot[['SFARI_070825', 'Fu_NDD664', 'Satterstrom_ASD102', 'Wang_NDD615']].apply(combine_or_no, axis=1)
        anno_reg_annot['Location'] = '.'
        anno_reg_annot.loc[anno_reg_annot['3UTR_genes'] != '', 'Location'] = '3UTR'
        anno_reg_annot.loc[anno_reg_annot['5UTR_genes'] != '', 'Location'] = '5UTR'
        anno_reg_annot.loc[anno_reg_annot['CDS_genes']  != '', 'Location'] = 'CDS'
        anno_reg_annot['Location'] = np.where((anno_reg_annot['Location'] == '.') & (anno_reg_annot['ENCODE_cCRE'].notna()), 'REG', anno_reg_annot['Location'])
        anno_reg_annot.to_csv(output.tsv, sep='\t', index=True)


