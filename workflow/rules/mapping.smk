## Map PE reads
rule bwa_mem:
    input:
        reads=get_fastq,
        idx=multiext("ref/NC_000962.3.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        bam="mapped/{sample}.bam",
        index="mapped/{sample}.bam.bai",
    log:
        "logs/bwa_mem_sambamba/{sample}.log",
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
    threads: config["BWA"]["threads"]
    wrapper:
        "v1.7.1/bio/bwa/mem-samblaster"
