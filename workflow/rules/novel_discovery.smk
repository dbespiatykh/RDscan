## Find regions with low coverage
rule compute_coverage:
    input:
        bam="mapped/{sample}.bam",
        mask=config["files"]["is6110"],
    output:
        bed="bed/{sample}.bed",
    conda:
        "../envs/calculations.yaml"
    log:
        "logs/bedtools/{sample}.coverage.log"
    shell:
        """
        (parallel bedtools genomecov -bga -ibam ::: {input.bam}\
         | awk '$4<'$(samtools depth  -aa {input.bam} |  awk '{{sum+=$3}} END {{print sum/4411532*0.1}}')''\
            | bedtools merge -d 1500 | bedtools subtract -f 0.30 -a stdin -b {input.mask} -A > {output.bed}) 2> {log}
        """


## Convert BED files to VCF
rule make_vcf:
    input:
        bed_cov=rules.compute_coverage.output.bed,
    output:
        vcf=temp("vcf/{sample}.raw.vcf"),
    conda:
        "../envs/calculations.yaml"
    log:
        "logs/survivor/{sample}.bedtovcf.log",
    shell:
        """
        SURVIVOR bedtovcf {input.bed_cov} DEL {output.vcf}
        """


## Remove unnecessary info from VCFs
rule clean_vcf:
    input:
        "vcf/{sample}.raw.vcf",
    output:
        temp("vcf/{sample}.vcf"),
    conda:
        "../envs/calculations.yaml"
    log:
        "logs/shell/{sample}.clean_vcf.log"
    shell:
        "(sed -e 's/;CIPOS=0,0;CIEND=0,0//' -e 's/-1/1/' -e 's/Sample/{wildcards.sample}/' {input} > {output}) 2> {log}"


## Compress VCFs
rule bgzip:
    input:
        "vcf/{sample}.vcf",
    output:
        "vcf/{sample}.vcf.gz",
    threads: 1
    log:
        "logs/bgzip/{sample}.log",
    wrapper:
        "v1.7.1/bio/bgzip"


## Index VCFs
rule bcftools_index:
    input:
        "vcf/{sample}.vcf.gz",
    output:
        "vcf/{sample}.vcf.gz.csi",
    log:
        "logs/bcftools/{sample}.index.log"
    wrapper:
        "v1.7.1/bio/bcftools/index"


rule filter_vcf:
    input:
        "vcf/{sample}.vcf.gz",
        "vcf/{sample}.vcf.gz.csi",
    output:
        "vcf/{sample}.filtered.vcf",
    log:
        "logs/bcftools/{sample}.filter.vcf.log",
    params:
        filter=config["filters"]["main"],
    wrapper:
        "v1.7.1/bio/bcftools/filter"


## Create list of VCF files
rule make_list:
    input:
        expand("vcf/{sample}.filtered.vcf", sample=samples.index),
    output:
        temp("vcf/vcf_list.txt"),
    conda:
        "../envs/calculations.yaml"
    log:
        "logs/shell/make_list.log",
    shell:
        "(ls {input} > {output}) 2> {log}"


## Merge VCF files
rule merge_vcf:
    input:
        "vcf/vcf_list.txt",
    output:
        "vcf/merged.vcf",
    conda:
        "../envs/calculations.yaml"
    log:
        "logs/survivor/merge.log",
    shell:
        "SURVIVOR merge {input} 5000 1 1 0 0 100 {output} 2> {log}"


## Convert multisample VCF to tab-separated format
rule convert2table:
    input:
        "vcf/merged.vcf",
    output:
        temp("RD_table.raw.tsv"),
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/gatk/variants2table.log",
    shell:
        "gatk VariantsToTable -V {input} -F CHROM -F POS -F END -F SVLEN -F SVTYPE -F LOF -GF LN -O {output} 2> {log}"


## Annotate putative RD regions
rule annotate_rds:
    input:
        "RD_table.raw.tsv",
        config["files"]["rds"],
    output:
        "results/RD_putative.tsv",
    log:
        "logs/annotate/novel.log",
    conda:
        "../envs/calculations.yaml"
    script:
        "../scripts/annotate_novel.py"
