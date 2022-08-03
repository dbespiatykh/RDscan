# Compute median coverage in known RD regions
rule mosdepth_bed:
    input:
        bam="mapped/{sample}.bam",
        bai="mapped/{sample}.bam.bai",
        bed=config["files"]["rds"],
    output:
        "mosdepth_bed/{sample}.mosdepth.global.dist.txt",
        "mosdepth_bed/{sample}.mosdepth.region.dist.txt",
        "mosdepth_bed/{sample}.regions.bed.gz",
        summary="mosdepth_bed/{sample}.mosdepth.summary.txt",
    log:
        "logs/mosdepth_bed/{sample}.log",
    params:
        extra="--no-per-base --use-median --fast-mode",
    threads: 4
    wrapper:
        "v1.7.1/bio/mosdepth"


## Calculate the proportion of the read depth in RD regions to the total chromosome depth
rule calculate_proportion:
    input:
        "mosdepth_bed/{sample}.regions.bed.gz",
    output:
        "bed/{sample}.proportion.bed",
    log:
        "logs/annotate/{sample}.proportion.log",
    conda:
        "../envs/calculations.yaml"
    script:
        "../scripts/proportions.py"


## Concatenate BED files
rule concatenate:
    input:
        bed=expand("bed/{sample}.proportion.bed", sample=samples.index),
    output:
        "results/all_concat.bed",
    log:
        "logs/annotate/concatenate.log",
    conda:
        "../envs/calculations.yaml"
    script:
        "../scripts/concatenate_bed.py"


## Generate tables with final RD calls
rule make_tables:
    input:
        "results/all_concat.bed",
    output:
        "results/RD_known.xlsx",
        "results/RD_known.tsv",
        "results/RD_known.bin.tsv",
    params:
        threshold=config["filters"]["threshold"],
    conda:
        "../envs/renv.yaml"
    log:
        "logs/annotate/r_tables.log",
    script:
        "../scripts/makeTables.R"
