## Map PE reads
rule bwa_mem:
    input:
        reads=get_fastq,
        idx=multiext(
            "resources/ref/NC_000962.3.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"
        ),
    output:
        bam="results/mapped/{sample}.bam",
        index="results/mapped/{sample}.bam.bai",
    log:
        "logs/bwa_mem_sambamba/{sample}.log",
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
    threads: config["BWA"]["threads"]
    wrapper:
        "v1.7.1/bio/bwa/mem-samblaster"
