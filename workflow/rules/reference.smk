rule download_genome:
    output:
        "ref/NC_000962.3.fa",
    params:
        accession=config["NCBI"]["H37Rv-reference-genome"],
    log:
        "logs/ncbi/ref_dwn.log",
    conda:
        "../envs/entrez.yaml"
    resources:
        ncbi_api_requests=1,
    shell:
        "((esearch -db nucleotide -query '{params.accession}' | "
        "efetch -format fasta > {output}) && [ -s {output} ]) 2> {log}"


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
