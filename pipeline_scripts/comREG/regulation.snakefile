#refer to: https://github.com/EichlerLab/svpop/blob/main/rules/variant/anno/regions.snakefile from SVpop.
import pandas as pd

sv_file = 'annotation_reg_input.bed'

manifest_df = pd.read_csv('manifest.txt', sep='\t', header=0, index_col="peak")
wildcard_constraints:
    peak='|'.join(manifest_df.index)

def get_peak(wildcards):
    return manifest_df.at[wildcards.peak, "path"]


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

rule variant_anno_max_constraint:
    input:
        bed=sv_file,
        nc_constraint='data/Non-coding_constraint/Supplementary_Data_2.bed'
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
        segdup='data/UCSC_repeat/SegDup.bed'
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


rule variant_anno_max_DHS_MEAN:
    input:
        bed=sv_file,
        dhs='data/dhs2020/dhs_200.bed'
    output:
        dhs_tsv='temp/dhs_sv.txt'
    resources:
        mem=10,
        hrs=24,
    envmodules:
        "miniconda/4.12.0"
    threads: 1
    #benchmark: "benchmark/oreganno.log"
    shell:
        """
        bedtools intersect -a {input.bed} -b {input.dhs} -wo|sort -k4,4 -k8,8r|awk '!seen[$4]++' |cut -f4,8 > {output.dhs_tsv}
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
        platinumTRs='temp/platinumTRs_sv.txt',
        RepeatMasker='temp/RepeatMasker_sv.txt',
        SimpleRepeat='temp/SimpleRepeat_sv.txt',
        tRNA='temp/tRNA_sv.txt',
        snRNA='temp/snRNA_sv.txt',
        lincRNA='temp/lincRNA_sv.txt',
        sno_miRNA='temp/sno_miRNA_sv.txt',
        dhs=rules.variant_anno_max_DHS_MEAN.output.dhs_tsv,
        segdup=rules.variant_anno_max_segdup.output.segdup_tsv,
        nc_constraint=rules.variant_anno_max_constraint.output.nc_constraint_tsv,
    output:
        tsv='allreg_sv.txt',
    resources:
        mem=10,
        hrs=24,
    threads: 1
    #benchmark: "benchmark/all.log"
    run:
        # Read SVs
        df = pd.read_csv(input.bed, sep='\t', squeeze=False)
        annotation_columns = ['ENCODE_cCREs', 'ENCODE_H3K27Ac', 'ENCODE_H3K4Me1', 'ENCODE_H3K4Me3', 'ENCODE_TFcluster', 'DHS', 'OREGANNO', 'geneHancer', 'CorticalMap', 'CUT&Tag_H3K27ac', 'CUT&Tag_H3K27me3', 'CUT&Tag_H3K4me3', 'CUT&Tag_CTCF', 'SegDup', 'RepeatMasker', 'SimpleRepeat', 'platinumTRs', 'ncConstraint', 'tRNA', 'snRNA', 'lincRNA', 'sno_miRNA']
        for col in annotation_columns:
            df[col] = None  # Initialize all columns with None

        df.set_index('ID', inplace=True)
        annotation_files = {
            'ENCODE_cCREs':input.ccre,
            'geneHancer':input.genehancer,
            'OREGANNO':input.oreganno,
            'CorticalMap':input.atac,
            'ENCODE_H3K27Ac':input.H3K27Ac,
            'ENCODE_H3K4Me1':input.H3K4Me1,
            'ENCODE_H3K4Me3':input.H3K4Me3,
            'ENCODE_TFcluster':input.TF,
            'CUT&Tag_H3K27ac':input.K_H3K27ac,
            'CUT&Tag_H3K27me3':input.K_H3K27me3,
            'CUT&Tag_H3K4me3':input.K_H3K4me3,
            'CUT&Tag_CTCF':input.K_CTCF,
            'platinumTRs':input.platinumTRs,
            'RepeatMasker':input.RepeatMasker,
            'SimpleRepeat':input.SimpleRepeat,
            'tRNA':input.tRNA,
            'snRNA':input.snRNA,
            'lincRNA':input.lincRNA,
            'sno_miRNA':input.sno_miRNA,
            'DHS':input.dhs,
            'SegDup':input.segdup,
            'ncConstraint':input.nc_constraint
            }

        for column, file_path in annotation_files.items():
            annotation_df = pd.read_csv(file_path, sep='\t', header=None, names=['ID', 'value'])
            annotation_df.set_index('ID', inplace=True)
            df.loc[df.index.isin(annotation_df.index), column] = annotation_df.loc[df.index.intersection(annotation_df.index), 'value']

        df.to_csv(output.tsv, sep='\t', index=True)

