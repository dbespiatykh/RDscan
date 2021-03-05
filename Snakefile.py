import os
import glob
import shutil
import datetime
import pandas as pd

os.makedirs("BAM", exist_ok=True)
os.makedirs("BED", exist_ok=True)
os.makedirs("VCF", exist_ok=True)

## ------------------------------------------------------------------------------------ ##
## Globals
## ------------------------------------------------------------------------------------ ##
SAMPLES, = glob_wildcards("reads/{smp}_1.fastq.gz")

## ------------------------------------------------------------------------------------ ##
## Workflow rule
## ------------------------------------------------------------------------------------ ##
rule all:
    input:
        expand("BAM/{smp}.bam", smp = SAMPLES),
        expand("BAM/{smp}.bam.bai", smp = SAMPLES),
        expand("BED/{smp}.bed", smp = SAMPLES),
        expand("VCF/{smp}_raw.vcf", smp = SAMPLES),
        expand("VCF/{smp}_cleaned.vcf", smp = SAMPLES),
        expand("VCF/{smp}_final.vcf", smp = SAMPLES),
        expand("VCF/{smp}_final.vcf.gz", smp = SAMPLES),
        expand("VCF/{smp}_final.vcf.gz.csi", smp = SAMPLES),
        expand("VCF/{smp}_final.duphold.vcf", smp = SAMPLES),
        expand("VCF/{smp}_final.filtered.vcf", smp = SAMPLES),
        'VCF/vcf_list.txt',
        "VCF/merged_raw.vcf",
        "VCF/merged.vcf",
        "RD_table_raw.tsv",
        "RD_table_cleaned.tsv",
        "RD_table.tsv",
        expand("BED/{smp}{ext}", smp = SAMPLES, ext = [".mosdepth.global.dist.txt",
        ".mosdepth.region.dist.txt",
        ".mosdepth.summary.txt",
        ".regions.bed.gz",
        ".regions.bed.gz.csi",
        ".regions.bed",
        ".proportion.bed"]
        ),
        "all_concat.bed",
        "RD_known.tsv",
        "RD_known.xlsx"

## ------------------------------------------------------------------------------------ ##
## RDScan rules
## ------------------------------------------------------------------------------------ ##
## Map reads to the H37Rv
rule mapping:
    input: r1 = "reads/{smp}_1.fastq.gz",
           r2 = "reads/{smp}_2.fastq.gz",
           ref= "db/H37Rv.fna"
    output: bam = "BAM/{smp}.bam"
    message: """--- Mapping reads ---"""
    shell: """
        bwa mem  -Y -M -R '@RG\\tID:{wildcards.smp}\\tSM:{wildcards.smp}' -t {config[threads]} {input.ref} {input.r1} {input.r2}\
         | samtools sort -@ 2 -o {output.bam} -
  """

## Generate index files for BAM alignments
rule index:
    input: bam = "BAM/{smp}.bam"
    output: bai = "BAM/{smp}.bam.bai"
    message: """--- Indexing ---"""
    shell: """
        samtools index {input.bam}
  """

## Find regions with low coverage
rule compute_coverage:
    input: bam = "BAM/{smp}.bam",
           mask = "db/IS6110.bed"
    output: bed = "BED/{smp}.bed"
    message: """--- Computing coverage ---"""
    shell: """
        parallel bedtools genomecov -bga -ibam ::: {input.bam}\
         | awk '$4<'$(samtools depth  -aa {input.bam} |  awk '{{sum+=$3}} END {{print sum/4411532*0.1}}')''\
            | bedtools merge -d 1500 | bedtools subtract -a stdin -b {input.mask} -A > {output.bed}
  """

## Convert BED files to VCF
rule make_vcf:
    input: bed_cov = rules.compute_coverage.output.bed
    output: vcf = "VCF/{smp}_raw.vcf"
    message: """--- Making VCF files ---"""
    shell: """
        SURVIVOR bedtovcf {input.bed_cov} DEL {output.vcf}
  """

## Remove unnecessary info from VCFs
rule clean_vcf:
    input: "VCF/{smp}_raw.vcf"
    output: "VCF/{smp}_cleaned.vcf"
    run:
      with open(input[0], 'r') as file :
        filedata = file.read()
        filedata = filedata.replace(';CIPOS=0,0;CIEND=0,0', '')
        filedata = filedata.replace('-1', '1')
      with open(output[0], 'w') as file:
        file.write(filedata)

## Add samplenames to VCFs
rule change_samplenames:
    input: "VCF/{smp}_cleaned.vcf"
    output: "VCF/{smp}_final.vcf"
    message: """--- Adding samplenames to VCF ---"""
    shell: """
        gatk RenameSampleInVcf -I {input} --NEW_SAMPLE_NAME {wildcards.smp} -O {output}
  """

## Compress VCF files
rule compress_vcf:
    input: "VCF/{smp}_final.vcf"
    output: "VCF/{smp}_final.vcf.gz"
    message: """--- Compressing VCF ---"""
    shell: """
        bcftools convert -Oz -o {output} {input}
  """

## Generate index files for VCFs
rule index_vcf:
    input: "VCF/{smp}_final.vcf.gz"
    output: "VCF/{smp}_final.vcf.gz.csi"
    message: """--- Indexing VCF ---"""
    shell: """
        bcftools index -f {input}
  """

## Calculate the depth of deletions flanking regions
rule estimate_depth:
    input: vcf = "VCF/{smp}_final.vcf.gz",
           csi = "VCF/{smp}_final.vcf.gz.csi",
           fasta = "db/H37Rv.fna",
           bam = "BAM/{smp}.bam",
           bai = "BAM/{smp}.bam.bai"
    output: "VCF/{smp}_final.duphold.vcf"
    conda: "envs/duphold.yaml"
    message: """--- Estimating depth ---"""
    shell: """
        duphold -v {input.vcf} -b {input.bam} -f {input.fasta} -o {output}
  """

## Filter VCF files
rule filter_vcf:
    input: "VCF/{smp}_final.duphold.vcf"
    output: "VCF/{smp}_final.filtered.vcf"
    message: """--- Filtering VCF ---"""
    shell: """
        bcftools view -i '(SVTYPE = "DEL" & FMT/DHFFC[0] < {config[DHFFC]} & INFO/SVLEN[0] > {config[minSVLEN]} & INFO/SVLEN[0] < {config[maxSVLEN]})' {input} > {output}
  """

## Create list of VCF files
rule make_list:
    input: expand("VCF/{smp}_final.filtered.vcf", smp = SAMPLES)
    output: "VCF/vcf_list.txt"
    message: """--- Making list of vcf files ---"""
    shell: 'ls {input} > {output}'

## Merge VCF files
rule merge_vcf:
    input: 'VCF/vcf_list.txt'
    output: "VCF/merged_raw.vcf"
    message: """--- Merging VCF files ---"""
    shell: """
        SURVIVOR merge {input} 5000 1 1 0 0 100 {output}
  """

## Annotate VCF files
rule annotate_vcf:
    input: 'VCF/merged_raw.vcf'
    output: "VCF/merged.vcf"
    message: """--- Annotating VCF ---"""
    shell: """
        snpEff ann m_tuberculosis_H37Rv -noLog -noStats -no-downstream -no-upstream -no-utr -download -o vcf {input} > {output}
  """

## Convert multisample VCF to tab-separated format
rule convert2table:
    input: "VCF/merged.vcf"
    output: "RD_table_raw.tsv"
    message: """--- Converting VCF ---"""
    shell: """
        gatk VariantsToTable -V {input} -F CHROM -F POS -F END -F SVLEN -F SVTYPE -F LOF -GF LN -O {output}
  """

## Remove unnecessary info from the table
rule clean_table:
    input: "RD_table_raw.tsv"
    output: "RD_table_cleaned.tsv"
    run:
         df = pd.read_csv(input[0], sep='\t', dtype='unicode')
         df['SVLEN'] = df['SVLEN'].str.replace(r'-', '')
         df.rename(columns={'POS': 'START', 'SVLEN': 'LENGTH', 'SVTYPE': 'TYPE'}, inplace=True)
         df.columns = df.columns.str.rstrip('.LN')
         df = df.replace({'0':''})
         df = df[df.CHROM != 'RvD1']
         df = df[df.CHROM != 'TbD1']
         lof = df['LOF']
         df.drop(labels=['LOF'], axis=1, inplace=True)
         lof = lof.str.extractall('(Rv\w+)').groupby(level=0)[0].apply(' '.join)
         df.insert(4, 'LOF', lof)
         df['LOF'] = df['LOF'].str.replace(r'\b(\w+)(\s+\1)+\b', r'\1', regex=True)
         df['LOF'] = df['LOF'].str.replace(" ",",")
         df.to_csv(output[0], sep='\t', index=False)

## Annotate putative RD regions
rule annotate_RD:
    input: "RD_table_cleaned.tsv",
           "db/RD_df.bed"
    output: "RD_table.tsv"
    run:
         df1 = pd.read_csv(input[0], sep='\t')
         df2 = pd.read_csv(input[1], sep='\t')
         s1 = df1.START.values
         e1 = df1.END.values
         s2 = df2.START.values
         e2 = df2.END.values
         r = pd.DataFrame(((e2-s1[:,None])>0)&((s2-e1[:,None])<0)).dot(df2.RD+'|').str[:-1]
         df1['RD'] = r.values
         rd = df1['RD']
         df1.drop(labels=['RD'], axis = 1, inplace=True)
         df1.insert(5, 'RD', rd)
         df1.to_csv('RD_table.tsv', sep='\t', index=False)

## Compute median coverage in known RD regions
rule compute_RD_coverage:
    input: bam = "BAM/{smp}.bam",
           bai = "BAM/{smp}.bam.bai",
           rd = "db/RD.bed"
    output: glo = "BED/{smp}.mosdepth.global.dist.txt",
            reg = "BED/{smp}.mosdepth.region.dist.txt",
            smr = "BED/{smp}.mosdepth.summary.txt",
            bedgz = "BED/{smp}.regions.bed.gz",
            ind = "BED/{smp}.regions.bed.gz.csi"
    params: prefix = "BED/{smp}"
    message: """--- Computing RD coverage ---"""
    shell: """
        parallel mosdepth -n -m -x -b {input.rd} {params.prefix} ::: {input.bam}
  """

## Decompress BED files
rule decompress:
    input: bedgz = "BED/{smp}.regions.bed.gz"
    output: bed = "BED/{smp}.regions.bed"
    threads: 2
    message: """--- Decompressing ---"""
    shell: """
        pigz -p {threads} -dc {input.bedgz} > {output.bed}
  """

## Calculate the proportion of the read depth in RD regions to the total chromosome depth
rule calculate_proportion:
    input: "BED/{smp}.regions.bed"
    output: "BED/{smp}.proportion.bed"
    run:
         df = pd.read_csv(input[0], sep='\t',
              names=['chrom', 'start', 'end', 'RD', 'depth'],
              dtype={'chrom':str, 'start':str, 'end':str, 'RD':str, 'depth':int})
         df.loc[:,'depth'] = df.loc[:,'depth'].div(df.iloc[0]['depth'])
         df['Sample'] = os.path.basename(input[0]).split('.regions.bed')[0]
         df = df[df.RD != 'H37Rv']
         df.to_csv(output[0], sep='\t', index=False)

## Concatenate BED files
rule concatenate:
    input: expand("BED/{smp}.proportion.bed", smp = SAMPLES)
    output: "all_concat.bed"
    run:
        files = glob.glob("BED/{smp}*.proportion.bed")
        df = pd.concat([pd.read_csv(fp, sep='\t') for fp in input])
        df = df.drop(['chrom', 'start', 'end'],  axis=1)
        df.to_csv(output[0], sep='\t', index=False)

## Generate tables with final RD calls
rule make_tables:
    input: "all_concat.bed"
    output: tab = "RD_known.tsv",
            xlsx = "RD_known.xlsx",
            bin = "RD_known.bin.tsv"
    message: """--- Making tables ---"""
    shell: "Rscript scripts/makeTables.R -T {config[threshold]} -i {input} -x {output.xlsx} -t {output.tab} -b {output.bin}"

onsuccess:
    dateStr = datetime.datetime.now().strftime("%m-%d-%Y_%Hh%Mm%Ss")
    resultsDir = "Results" + "_" + dateStr
    os.makedirs(resultsDir, exist_ok=True)
    shutil.move("RD_known.bin.tsv", resultsDir)
    shutil.move("RD_known.tsv", resultsDir)
    shutil.move("RD_known.xlsx", resultsDir)
    shutil.move("RD_table.tsv", resultsDir)
    shutil.rmtree(".snakemake", ignore_errors=True)
    shutil.rmtree("VCF", ignore_errors=True)
    shutil.rmtree("BED", ignore_errors=True)
    os.remove("RD_table_raw.tsv")
    os.remove("RD_table_cleaned.tsv")
    os.remove("all_concat.bed")
    print("Completed successfully")

onerror:
    print("An error occurred")