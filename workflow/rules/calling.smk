## Call variants
rule haplotype_caller:
    input:
        bam="BAM/{smp}.bam",
        idx="BAM/{smp}.bam.bai",
        ref="ref/NC_000962.3.fa",
        fai="ref/NC_000962.3.fa.fai",
        dic="ref/NC_000962.3.dict",
    output:
        gvcf="VCF/{smp}.g.vcf.gz",
    log:
        "logs/gatk/haplotypecaller/{smp}.log",
    params:
        extra="-ploidy 1 -mbq 20",
    threads: config["hapcall_threads"]
    resources:
        mem_mb=config["hapcall_mem"],
    wrapper:
        "v1.7.1/bio/gatk/haplotypecaller"


## Import VCFs to GenomicsDB
rule genomics_db_import:
    input:
        gvcfs=expand("VCF/{smp}.g.vcf.gz", smp=SAMPLES),
    output:
        db=directory("db"),
    log:
        "logs/gatk/genomicsdbimport.log",
    params:
        intervals="NC_000962.3",
        db_action="create",
        extra="--batch-size 100",
    resources:
        mem_mb=config["gendb_mem"],
    wrapper:
        "v1.7.1/bio/gatk/genomicsdbimport"


## Genortype Variants
rule genotype_gvcfs:
    input:
        genomicsdb="db",
        ref="ref/NC_000962.3.fa",
    output:
        vcf="VCF/all.vcf.gz",
    log:
        "logs/gatk/genotypegvcfs.log",
    params:
        extra="-ploidy 1",
    resources:
        mem_mb=config["genotype_mem"],
    wrapper:
        "v1.7.1/bio/gatk/genotypegvcfs"


## Filter SNPs
rule gatk_filter:
    input:
        vcf="VCF/all.vcf.gz",
        ref="ref/NC_000962.3.fa",
    output:
        vcf="VCF/all.filtered.vcf.gz",
    log:
        "logs/gatk/filter/snps.filter.log",
    params:
        filters={"Failfilter": "QD < 2.0 || DP < 10 || FS > 60.0 || MQ < 40.0"},
    resources:
        mem_mb=config["filter_mem"],
    wrapper:
        "v1.7.1/bio/gatk/variantfiltration"


## Select only PASS variants
rule gatk_select_passed:
    input:
        vcf="VCF/all.filtered.vcf.gz",
        ref="ref/NC_000962.3.fa",
    output:
        vcf="VCF/all.snps.pass.vcf.gz",
    log:
        "logs/gatk/select/snps.pass.log",
    params:
        extra="--exclude-filtered --select-type-to-include SNP",
    resources:
        mem_mb=config["select_mem"],
    wrapper:
        "v1.7.1/bio/gatk/selectvariants"


## Variants to table
rule gatk_variants_to_table:
    input:
        vcf="VCF/all.snps.pass.vcf.gz",
    output:
        tab="Results/all.snps.pass.txt",
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/gatk/vartotable/vars2table.log",
    shell:
        "gatk VariantsToTable \
        -V {input.vcf} \
        -O {output.tab} \
        -F POS \
        -F REF \
        -GF GT 2> {log}"
