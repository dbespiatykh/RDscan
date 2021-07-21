<img align ="left" src=img/RDscan_logo.png width=250px style="padding-right: 25px; padding-top: 25px;">

# pipeline for MTBC putative regions of difference discovery

[![citation](https://img.shields.io/badge/DOI-10.1128%2FmSphere.00535--21-9f1d21)](https://doi.org/10.1128/mSphere.00535-21)

- [Description](#Description)
- [Installation](#Installation)
- [Usage](#Usage)
- [Output](#Output)
- [Citation](#Citation)

## Description

RDscan is a [snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline to find deletions and putative [regions of difference](https://jb.asm.org/content/178/5/1274.short) (RD) in [mycobacterium tuberculosis complex](https://en.wikipedia.org/wiki/Mycobacterium_tuberculosis_complex) (MTBC) genomes, it is also capable to determine already known or user defined RDs.

## Installation

Use the [Conda](https://docs.conda.io/en/latest/) package manager and [BioConda](https://bioconda.github.io/index.html) channel to install RDscan.

If you do not have conda installed do the following:

```bash
# Download Conda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# Set permissions
chmod -X Miniconda3-latest-Linux-x86_64.sh
# Install
bash Miniconda3-latest-Linux-x86_64.sh
```

Set up channels:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Get **RDscan** pipeline:

```bash
git clone https://github.com/dbespiatykh/RDscan.git
```

Install all required dependencies:

```bash
cd RDscan
conda env create --file environment.yml
```

## Usage

#### Rulegraph of the pipeline

<br>

![Rulegraph](img/Rulegraph.png)
<br>
:point_right: In project folder make `reads` folder and move to it your paired-end [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) reads (suffix should be: `_1.fastq.gz` and `_2.fastq.gz`).

:grey_exclamation: You can also add your own RDs positions, which are not present in the `db/RD.bed` table.\
Columns are:\
**1.** Chromosome name;\
**2.** Start position of the RD;\
**3.** End position of the RD;\
**4.** Name of the RD.

|           |       |       |              |
| --------- | ----- | ----- | ------------ |
| NC_000962 | 29988 | 34322 | RDcap_Spain1 |
| NC_000962 | 34789 | 35209 | RD301        |
| NC_000962 | 76163 | 84826 | RD105ext     |
| NC_000962 | 79571 | 83036 | RD105        |
| …         | …     | …     | …            |

:file_folder: Project folder should have the following structure:

```bash
📂RDscan/
├── Snakefile
├── config.json
├── db
│   ├── H37Rv.fna
│   ├── H37Rv.fna.amb
│   ├── H37Rv.fna.ann
│   ├── H37Rv.fna.bwt
│   ├── H37Rv.fna.fai
│   ├── H37Rv.fna.pac
│   ├── H37Rv.fna.sa
│   ├── IS6110.bed
│   ├── RD.bed
│   └── RD_df.bed
├── envs
│   └── duphold.yaml
├── reads
│   ├── sample-1_1.fastq.gz
│   ├── sample-1_2.fastq.gz
│   ├── sample-2_1.fastq.gz
│   ├── sample-2_2.fastq.gz
│   └── ...
└── scripts
    └── makeTables.R

```

<br>

Activate **RDscan** environment:

```bash
conda activate RDscan
```

Run pipeline:

```bash
snakemake -j {Number of cores} --use-conda --configfile config.json
```

It is recommended to use dry run if you are running pipeline for the first time, to see if everything is in working order, for this you can use `-n` flag:

```bash
snakemake -j {Number of cores} --use-conda --configfile config.json -n
```

Parameters that can be adjusted in `config.json`:

```json
{
  "threads": 8,
  "threshold": 0.05,
  "DHFFC": 0.1,
  "minSVLEN": 200,
  "maxSVLEN": 30000
}
```

- `threads` - number of threads to use for [BWA-MEM](https://github.com/lh3/bwa) (8 by default);

> e.g. If you run `snakemake -j 16` and `threads` parameter in `config.json` is equals to `2` , then snakemake will execute up to `8` instances of the `BWA-MEM` in the mapping rule.

- `threshold` - threshold value to use for coverage condition filtering (0.05 by default)
- `DHFFC` - [Duphold](https://github.com/brentp/duphold) flank fold-change (0.1 by default);
- `minSVLEN` - minimum deletion length in bp (200 by default);
- `maxSVLEN` - maximum deletion length in bp (30,000 by default)

## Output

Output in the `Results_MM-DD-YYYY_HHh-MMm-SSs` directory will contain four tables: `RD_table.tsv`, `RD_known.tsv`, `RD_known.xlsx` and `RD_known.bin.tsv`

Example of the `RD_table.tsv`:
Table containing all discovered putative RDs.\
**RD** - Known RDs that intersects with deletion breakpoints;\
**LOF** - loss of function genes annotated with [snpEff](https://pcingola.github.io/SnpEff/);\
**LENGTH** - Estimated size of predicted deletion.
Values in cells represent deletion length in the sample.
| CHROM | START | END | LENGTH | LOF | RD | TYPE | ERR015582 | ERR017778 | ERR017782 | ERR019852 |
| ---------- | ------ | ------ | ------ | --------------- | --- | ---- | --------- | --------- | --------- | --------- |
| NC_000962 | 333828 | 338580 | 5800 | Rv0278c,Rv0279c | | DEL | 7113 | 7084 | 7050 |
| NC_000962 | 340400 | 340645 | 245 | Rv0280 | | DEL | | | | |
| NC_000962 | 350935 | 351175 | 238 | | | DEL | | 300 | 204 | 240 |
| NC_000962 | 361769 | 362988 | 1391 | Rv0297 | | DEL | 1833 | 1392 | 1833 | 1390 |

Example of the `RD_known.tsv`:
Table containing proportion of coverage in particular RDs.
| Sample | N-RD25_tbA | N-RD25_tbB | N-RD25bov/cap | N-RD25das |
| --------- | ----------- | ----------- | ------------- | --------- |
| ERR015582 | 0.883562 | 0.856164 | 0.856164 | 0.808219 |
| ERR017778 | 0 | 0 | 0 | 0.41791 |
| ERR017782 | 1.021277 | 1.042553 | 1.106383 | 0.978723 |
| ERR019852 | 0 | 0 | 0 | 0.386364 |

Example of the `RD_known.xlsx`:
Same as the `RD_known.tsv`, but in a [XLSX](https://en.wikipedia.org/wiki/Microsoft_Excel) format with applied contiditional formatting.\
Conditional formatting corresponds with threshold value in a `config.json` file.

![](img/RD_known.xlsx.png)

`RD_known.bin.tsv` is the same as `RD_known.tsv` and `RD_known.xlsx`, but in binary form:
| Sample | N-RD25_tbA | N-RD25_tbB | N-RD25bov/cap | N-RD25das |
| --------- | ----------- | ----------- | ------------- | --------- |
| ERR015582 | 0 | 0 | 0 | 0 |
| ERR017778 | 1 | 1 | 1 | 0 |
| ERR017782 | 0 | 0 | 0 | 0 |
| ERR019852 | 1 | 1 | 1 | 0 |

## Citation

If you use `RDscan` for your research, please cite the pipeline:

> D. Bespiatykh, J. Bespyatykh, I. Mokrousov, and E. Shitikov, A Comprehensive Map of Mycobacterium tuberculosis Complex Regions of Difference, mSphere, Volume 6, Issue 4, 21 July 2021, Page e00535-21, https://doi.org/10.1128/mSphere.00535-21

All references for the tools utilized by the `RDscan` can be found in the [`CITATIONS.md`](CITATIONS.md) file.

## License

[MIT](LICENSE)
