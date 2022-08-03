# General settings

Following the explanations in the `config.yml`, you should modify it according to your needs.

# Samples sheet

The location of this sheet must be specified in the `config.yml`.

It should be formatted like this:

| Run_accession | R1                             | R2                             |
| ------------- | ------------------------------ | ------------------------------ | --- |
| SRR2024996    | /path/to/SRR2024996_1.fastq.gz | /path/to/SRR2024996_2.fastq.gz |
| SRR2024925    | /path/to/SRR2024925_1.fastq.gz | /path/to/SRR2024925_2.fastq.gz |     |

**Run_accession** - Run accession number or sample name;\
**R1** - Path to the first read pair;\
**R2** - Path to the second read pair.

- Both `R1` and `R2`should be specified with reads paths.

# RDs

:point_right: You can add your own RD positions, which are not present in the `RD.bed` table.\
Columns are:

**1.** Chromosome name;\
**2.** Start position of the RD;\
**3.** End position of the RD;\
**4.** Name of the RD.

|             |       |       |              |
| ----------- | ----- | ----- | ------------ |
| NC_000962.3 | 29988 | 34322 | RDcap_Spain1 |
| NC_000962.3 | 34789 | 35209 | RD301        |
| NC_000962.3 | 76163 | 84826 | RD105ext     |
| NC_000962.3 | 79571 | 83036 | RD105        |
| …           | …     | …     | …            |
