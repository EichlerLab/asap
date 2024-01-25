import pandas as pd
import os

manifest_df = pd.read_table(
    config["manifest"], keep_default_na=False, dtype="str"
).set_index(
    ["sample"],drop=False
)  # sample,fofn


target_ref = config["target_ref"]
tech = config["tech"]

def calc_mem_gb(wildcards, input, attempt, threads):
    mb = max(1.5 * input.size_mb, 1000)
    gb = int(mb / 1000)

    if threads != 1:
        gb = int(max(gb / threads, 2))

    return gb * attempt

def get_reference(wildcards):
    reference_path = config["reference"][wildcards.ref]
    return reference_path


def get_y_estimates(wildcards):
    cell_names = get_fofn_df(sample_name=wildcards.sample).index.tolist()
    return expand(
        "results/{tech}/{ref}/{{sample}}/tmp/{{sample}}_{cell}-Ystats.tsv",
        tech=tech, ref=target_ref, cell=cell_names
    )


def get_fofn_df(sample_name):
    fp = manifest_df.at[sample_name, "file_path"]
    fofn_df = pd.read_table(fp,header=None,names=["file_path"])
    fofn_df["basename"] = fofn_df.apply(
        lambda row: os.path.basename(row.file_path),axis=1
    )
    fofn_df["n"] = fofn_df.index
    fofn_df["cell"] = fofn_df.basename.str.extract(r"(.*)\.fastq.*")

    # If no duplicates are present in file_path but exists in cell, then add unique identifier to cell
    if fofn_df["cell"].duplicated().any() and not fofn_df["file_path"].duplicated().any():
        fofn_df["cell"] = fofn_df.apply(lambda x: f"{x.cell}-{x.n}",axis=1)

    fofn_df.set_index(["cell"],inplace=True)
    return fofn_df

def get_cell_fastq(wildcards):
    fofn_df = get_fofn_df(sample_name=wildcards.sample)
    return fofn_df.at[wildcards.cell, "file_path"]

rule all:
    input:
        expand(
            "results/{tech}/{ref}/{sample}/{sample}-Ystats.tsv.gz",
            sample=manifest_df.index, tech=tech, ref=target_ref
        )

#####################
### HiFi-Specific ###
#####################
rule pbmm2_index:
    input:
        reference = get_reference
    output:
        reference_mmi="results/hifi/{ref}/resources/pbmm2/{ref}.mmi",
    wildcard_constraints:
        tech = "hifi"
    threads: 8
    resources:
        mem=calc_mem_gb,
        hrs=72,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "pbconda/202307",
    shell:
        """
        pbmm2 index --preset HiFi {input.reference} {output.reference_mmi}
        """

rule pbmm2:
    input:
        fastq = get_cell_fastq,
        reference = "results/hifi/{ref}/resources/pbmm2/{ref}.mmi"
    output:
        bam=temp("results/hifi/{ref}/{sample}/tmp/{sample}_{cell}.bam"),
        bai=temp("results/hifi/{ref}/{sample}/tmp/{sample}_{cell}.bam.bai"),
    threads: 8
    resources:
        mem=lambda wildcards, attempt: attempt * 4,
        hrs=72,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "pbconda/202307",
    shell:
        """
        pbmm2 align --preset HIFI {input.reference} {input.fastq} {output.bam} --sort --rg "@RG\tID:{wildcards.cell}\tSM:{wildcards.sample}\tPL:PACBIO"
        """

#####################
### ONT-Specific ###
#####################
rule align:
    input:
        fastq=get_cell_fastq,
        reference=get_reference,
    output:
        sam=temp("results/ont/{ref}/{sample}/tmp/{sample}_{cell}.sam"),
    params:
        mm2_params = MINIMAP2_PARAMS
    threads: 12
    resources:
        mem=lambda wildcards, attempt: attempt * 4,
        hrs=72,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "minimap2/2.26"
    shell:
        """
        minimap2 \
            -t {threads} \
            -L --secondary=no \
            -ax map-ont \
            {params.mm2_params} \
            {input.reference} {input.fastq} > {output.sam}
        """

rule index_align:
    input:
        sam="results/ont/{ref}/{sample}/tmp/{sample}_{cell}.sam",
    output:
        bam=temp("results/ont/{ref}/{sample}/tmp/{sample}_{cell}.bam"),
        bai=temp("results/ont/{ref}/{sample}/tmp/{sample}_{cell}.bam.bai"),
    threads: 12
    resources:
        mem=lambda wildcards, attempt: attempt * 2,
        hrs=72,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "samtools/1.14"
    shell:
        """
        samtools view -@ {threads} -b -o {output.bam} {input.sam} && samtools index -@ {threads} {output.bam}
        """

rule estimate_sex:
    input:
        bam = "results/{tech}/{ref}/{sample}/tmp/{sample}_{cell}.bam",
        bai = "results/{tech}/{ref}/{sample}/tmp/{sample}_{cell}.bam.bai"
    output:
        stats = temp("results/{tech}/{ref}/{sample}/tmp/{sample}_{cell}-stats.tsv"),
        ystats = temp("results/{tech}/{ref}/{sample}/tmp/{sample}_{cell}-Ystats.tsv")
    params:
        prefix = "results/{tech}/{ref}/{sample}/tmp/{sample}_{cell}"
    threads: 8
    resources:
        mem=lambda wildcards, attempt: attempt * 10,
        hrs=72,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "miniconda/4.12.0",
    shell:
        """
        ./sex-estimator.py {input.bam} --threads {threads} --prefix {params.prefix}
        """

rule merge_estimates:
    input:
        y_estimates = get_y_estimates
    output:
        merged = "results/{tech}/{ref}/{sample}/{sample}-Ystats.tsv.gz"
    threads: 1
    resources:
        mem=calc_mem_gb,
        hrs=72,
    run:
        df = pd.concat([pd.read_table(x, header=0, dtype=str) for x in input.y_estimates])
        df.to_csv(output.merged, header=True, index=False, sep="\t")
