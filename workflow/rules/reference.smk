from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider

NCBI = NCBIRemoteProvider(email=config["NCBI"]["email"])


rule download_genome:
    input:
        genome=NCBI.remote("NC_000962.3.fasta", db="nuccore"),
        rd=config["files"]["rvd1_tbd1"],
    output:
        "ref/NC_000962.3.fa",
    run:
        shell("cat {input.genome} {input.rd} > {output}")


rule samtools_genome_index:
    input:
        "ref/NC_000962.3.fa",
    output:
        "ref/NC_000962.3.fa.fai",
    log:
        "logs/samtools/ref_index.log",
    wrapper:
        "v1.7.1/bio/samtools/faidx"


rule bwa_index:
    input:
        "ref/NC_000962.3.fa",
    output:
        idx=multiext("ref/NC_000962.3.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa/ref_index.log",
    params:
        algorithm="bwtsw",
    wrapper:
        "v1.7.1/bio/bwa/index"
