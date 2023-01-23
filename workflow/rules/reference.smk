rule download_genome:
    output:
        "resources/ref/NC_000962.3.fa",
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
        "resources/ref/NC_000962.3.fa",
    output:
        "resources/ref/NC_000962.3.fa.fai",
    log:
        "logs/samtools/ref_index.log",
    wrapper:
        "v1.21.4/bio/samtools/faidx"


rule bwa_index:
    input:
        "resources/ref/NC_000962.3.fa",
    output:
        idx=multiext(
            "resources/ref/NC_000962.3.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"
        ),
    log:
        "logs/bwa/ref_index.log",
    params:
        algorithm="bwtsw",
    wrapper:
        "v1.21.4/bio/bwa/index"
