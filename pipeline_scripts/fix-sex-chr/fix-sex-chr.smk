import pandas as pd


configfile: "config/config.yaml"


MANIFEST = config.get("manifest", "config/manifest.tab")
MINIMAP_PARAMS = config.get("minimap_params", "-x asm20 --secondary=no -s 25000")
REF = config.get("reference")
MATCH_ARGS = config.get("match_args", "-f 0.95 -r")
COMP_ARGS = config.get("complement_args", "-f 0.95 -r -v")

manifest_df = pd.read_csv(MANIFEST, sep="\t", index_col="sample", comment="#")
haps = manifest_df.columns
sex_chr = ["chrY", "chrX"]


wildcard_constraints:
    sample="|".join(manifest_df.index),
    hap="|".join(manifest_df.columns),
    chr="|".join(sex_chr),


def get_hap_fasta(wildcards):
    return manifest_df.at[wildcards.sample, wildcards.hap]


def get_both_subset_paf(wildcards):
    paf_subset = {}
    for hap in haps:
        search_path = f"results/{wildcards.sample}/{wildcards.sample}_{hap}-subset_{wildcards.chr}.bed"
        paf_subset[hap] = search_path
    return paf_subset


def get_mend_fasta_inputs(wildcards):
    all_inputs = {}
    s = wildcards.sample
    # Meant to return matchs and complements for each chromosome
    # as well as headers and fasta for each haplotype
    # In addition to the hap1-chrX bed and the hap2-chrY bed to check for overlap
    for c, h in zip(sex_chr, haps):
        all_inputs[f"matches-{c}"] = f"results/{s}/{c}-matches.txt"
        all_inputs[f"complements-{c}"] = f"results/{s}/{c}-complements.txt"
        all_inputs[f"headers-{h}"] = f"results/{s}/{s}_{h}-fasta_headers.txt"
        all_inputs[f"fasta-{h}"] = manifest_df.at[s, h]
        all_inputs[f"bed-{h}-{c}"] = f"results/{s}/{s}_{h}-subset_{c}.bed"
    return all_inputs


rule all:
    input:
        expand(
            "results/{sample}/{sample}_{hap}.{ext}",
            sample=manifest_df.index,
            hap=haps,
            ext=["fasta", "fasta.fai"],
        ),


rule reassign_tab:
    input:
        expand(
            "results/{sample}/{sample}_reassign.tab",
            sample=manifest_df.index,
        ),


rule make_paf:
    input:
        fa=get_hap_fasta,
        ref=REF,
    output:
        paf="results/{sample}/{sample}_{hap}.paf",
    singularity:
        "docker://eichlerlab/align-basics:0.1"
    threads: 8
    resources:
        mem=12,
        hrs=72,
    shell:
        """
        minimap2 -c -t {threads} -K {resources.mem}G --eqx --cs {MINIMAP_PARAMS} {input.ref} {input.fa} -o {output.paf}
        """


rule extract_fasta_headers:
    input:
        fa=get_hap_fasta,
    output:
        headers=temp("results/{sample}/{sample}_{hap}-fasta_headers.txt"),
    threads: 1
    resources:
        mem=4,
        hrs=8,
    shell:
        """
        if ( file $( readlink -f {input.fa} ) | grep -q zip ); then
            zcat {input.fa} | grep '>' | cut -f2 -d'>' > {output.headers}
        else
            grep '>' {input.fa} | cut -f2 -d'>' > {output.headers}
        fi
        """


rule subset_paf:
    input:
        paf_file=rules.make_paf.output.paf,
    output:
        paf_bed=temp("results/{sample}/{sample}_{hap}-subset_{chr}.bed"),
    threads: 1
    resources:
        mem=4,
        hrs=12,
    shell:
        """
        cat {input.paf_file} | grep {wildcards.chr} | awk -F'\\t' '{{print $6,$8,$9,$1,$3,$4}}' OFS='\\t' > {output.paf_bed} || touch {output.paf_bed}
        """


rule extract_relevant_contigs:
    input:
        unpack(get_both_subset_paf),
    output:
        matches=temp("results/{sample}/{chr}-matches.txt"),
        complements=temp("results/{sample}/{chr}-complements.txt"),
    singularity:
        "docker://eichlerlab/binf-basics:0.1"
    threads: 1
    resources:
        mem=4,
        hrs=12,
    params:
        match_arg=MATCH_ARGS,
        comp_arg=COMP_ARGS,
    shell:
        """
        if [[ {wildcards.chr} == "chrX" ]]; then
            bedtools intersect \
                {params.match_arg} \
                -a {input.hap1} \
                -b {input.hap2} \
                | cut -f4 | sort -u > {output.matches}

            bedtools intersect \
                {params.comp_arg} \
                -a {input.hap1} \
                -b {input.hap2} \
                | cut -f4 | sort -u > {output.complements}
        else
            bedtools intersect \
                {params.match_arg} \
                -b {input.hap1} \
                -a {input.hap2} \
                | cut -f4 | sort -u > {output.matches}

            bedtools intersect \
                {params.comp_arg} \
                -b {input.hap1} \
                -a {input.hap2} \
                | cut -f4 | sort -u > {output.complements}
        fi
        """


rule mend_fasta:
    input:
        unpack(get_mend_fasta_inputs),
    output:
        chrX_add_to_hap2=temp("results/{sample}/chrX-add_to_hap2.fasta"),
        chrX_new_hap1=temp("results/{sample}/chrX-new_hap1.fasta"),
        chrY_add_to_hap1=temp("results/{sample}/chrY-add_to_hap1.fasta"),
        chrY_new_hap2=temp("results/{sample}/chrY-new_hap2.fasta"),
        reassign_tab="results/{sample}/{sample}_reassign.tab",
    singularity:
        "docker://eichlerlab/binf-basics:0.1"
    threads: 1
    resources:
        mem=4,
        hrs=24,
    shell:
        """
        # we want to move all chrY relevant contigs to hap1
        # we want to move all chrX relevant contigs to hap2
        ## Unless there is no overlap
        echo -e "ctg\\treason\\tchrom\\tdest" > {output.reassign_tab}
        if [[ $( cat {input.bed-hap1-chrY} | wc -l ) == 0 ]] && [[ $( cat {input.bed-hap2-chrX} | wc -l ) == 0 ]]; then
            touch {output.chrX_add_to_hap2}
            touch {output.chrY_add_to_hap1}
            ln -s $( readlink -f {input.fasta-hap1} ) {output.chrX_new_hap1}
            ln -s $( readlink -f {input.fasta-hap2} ) {output.chrY_new_hap2}
        else
            # extract chrX complements from hap1 fasta and then add to hap2 fasta
            seqtk subseq {input.fasta-hap1} {input.complements-chrX} | sed -E 's/(>)(.*)/\\1\\2-adjusted/g; s/h1/h2/g' | seqtk seq -l 80 > {output.chrX_add_to_hap2}
            sed 's/$/ complement chrX hap2/g' {input.complements-chrX} | awk -vOFS="\\t" '{{print $1,$2,$3,$4}}' 2>&1 | tee -a {output.reassign_tab}

            # remove the matches + complements from hap1 fasta
            seqtk subseq {input.fasta-hap1} <(grep -Ev "$(cat {input.matches-chrX} {input.complements-chrX})" {input.headers-hap1}) | seqtk seq -l 80 >  {output.chrX_new_hap1}
            sed 's/$/ duplicate chrX hap2/g' {input.matches-chrX} | awk -vOFS="\\t" '{{print $1,$2,$3,$4}}' 2>&1 | tee -a {output.reassign_tab}

            # extract chrY complements from hap2 fasta and then add to hap1 fasta
            seqtk subseq {input.fasta-hap2} {input.complements-chrY} | sed -E 's/(>)(.*)/\\1\\2-adjusted/g; s/h2/h1/g' | seqtk seq -l 80 > {output.chrY_add_to_hap1}
            sed 's/$/ complement chrY hap1/g' {input.complements-chrY} | awk -vOFS="\\t" '{{print $1,$2,$3,$4}}' 2>&1 | tee -a {output.reassign_tab}

            # remove the matches + complements from hap2 fasta
            seqtk subseq {input.fasta-hap2} <(grep -Ev "$(cat {input.matches-chrY} {input.complements-chrY})" {input.headers-hap2}) | seqtk seq -l 80 >  {output.chrY_new_hap2}    
            sed 's/$/ duplicate chrY hap1/g' {input.matches-chrY} | awk -vOFS="\\t" '{{print $1,$2,$3,$4}}' 2>&1 | tee -a {output.reassign_tab}
        fi 
        """


rule modified_fasta:
    input:
        chrX_add_to_hap2="results/{sample}/chrX-add_to_hap2.fasta",
        chrX_new_hap1="results/{sample}/chrX-new_hap1.fasta",
        chrY_add_to_hap1="results/{sample}/chrY-add_to_hap1.fasta",
        chrY_new_hap2="results/{sample}/chrY-new_hap2.fasta",
    output:
        hap1="results/{sample}/{sample}_hap1.fasta",
        hap2="results/{sample}/{sample}_hap2.fasta",
    threads: 1
    resources:
        mem=4,
        hrs=24,
    shell:
        """
        cat {input.chrX_new_hap1} {input.chrY_add_to_hap1} > {output.hap1} 
        cat {input.chrY_new_hap2} {input.chrX_add_to_hap2} > {output.hap2}
        """


rule faidx_fasta:
    input:
        hap1=rules.modified_fasta.output.hap1,
        hap2=rules.modified_fasta.output.hap2,
    output:
        hap1_fai="results/{sample}/{sample}_hap1.fasta.fai",
        hap2_fai="results/{sample}/{sample}_hap2.fasta.fai",
    threads: 1
    resources:
        mem=4,
        hrs=24,
    singularity:
        "docker://eichlerlab/binf-basics:0.1"
    shell:
        """
        samtools faidx {input.hap1}
        samtools faidx {input.hap2}
        """
