import os
import glob
import time
import shutil
import datetime
import pandas as pd

t_start = time.time()

os.makedirs("BAM", exist_ok=True)
os.makedirs("BED", exist_ok=True)
os.makedirs("VCF", exist_ok=True)
os.makedirs("logs", exist_ok=True)

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
        "RD_known.tsv",
        "RD_known.xlsx",
        "RD_known.bin.tsv",
        "RD_table.tsv"

## ------------------------------------------------------------------------------------ ##
## RDScan rules
## ------------------------------------------------------------------------------------ ##

## Map reads to the H37Rv
rule mapping:
    input: r1 = "reads/{smp}_1.fastq.gz",
           r2 = "reads/{smp}_2.fastq.gz",
           ref= "db/H37Rv.fna"
    output: bam = "BAM/{smp}.bam"
    log: "logs/bwa_{smp}.log"
    message: """--- Mapping {wildcards.smp} to H37Rv ---"""
    shell: """
        (bwa mem  -Y -M -R '@RG\\tID:{wildcards.smp}\\tSM:{wildcards.smp}' -t {config[threads]} {input.ref} {input.r1} {input.r2}\
         | samtools sort -@ 2 -o {output.bam} -) 2> {log}
  """

## Generate index files for BAM alignments
rule index:
    input: bam = "BAM/{smp}.bam"
    output: bai = "BAM/{smp}.bam.bai"
    message: """--- Indexing {wildcards.smp}.bam ---"""
    shell: """
        samtools index {input.bam}
  """

## Find regions with low coverage
rule compute_coverage:
    input: bam = "BAM/{smp}.bam",
           mask = "db/IS6110.bed"
    output: bed = "BED/{smp}.bed"
    message: """--- Computing coverage of {wildcards.smp}.bam ---"""
    shell: """
        parallel bedtools genomecov -bga -ibam ::: {input.bam}\
         | awk '$4<'$(samtools depth  -aa {input.bam} |  awk '{{sum+=$3}} END {{print sum/4411532*0.1}}')''\
            | bedtools merge -d 1500 | bedtools subtract -f 0.30 -a stdin -b {input.mask} -A > {output.bed}
  """

## Convert BED files to VCF
rule make_vcf:
    input: bed_cov = rules.compute_coverage.output.bed
    output: vcf = "VCF/{smp}_raw.vcf"
    message: """--- Making {wildcards.smp}_raw.vcf ---"""
    shell: """
        SURVIVOR bedtovcf {input.bed_cov} DEL {output.vcf}
  """

## Remove unnecessary info from VCFs
rule clean_vcf:
    input: "VCF/{smp}_raw.vcf"
    output: "VCF/{smp}_cleaned.vcf"
    message: """--- Cleaning {wildcards.smp}_raw.vcf ---"""
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
    log: "logs/gatk_addnames_{smp}.log"
    message: """--- Adding samplename to {wildcards.smp}_cleaned.vcf ---"""
    shell: """
        gatk RenameSampleInVcf -I {input} --NEW_SAMPLE_NAME {wildcards.smp} -O {output} 2> {log}
  """

## Compress VCF files
rule compress_vcf:
    input: "VCF/{smp}_final.vcf"
    output: "VCF/{smp}_final.vcf.gz"
    log: "logs/bcftools_comp_{smp}.log"
    message: """--- Compressing {wildcards.smp}.vcf ---"""
    shell: """
        bcftools convert -Oz -o {output} {input} 2> {log}
  """

## Generate index files for VCFs
rule index_vcf:
    input: "VCF/{smp}_final.vcf.gz"
    output: "VCF/{smp}_final.vcf.gz.csi"
    log: "logs/bcftools_index_{smp}.log"
    message: """--- Indexing VCF ---"""
    shell: """
        bcftools index -f {input} 2> {log}
  """

## Calculate the depth of deletions flanking regions
rule estimate_depth:
    input: vcf = "VCF/{smp}_final.vcf.gz",
           csi = "VCF/{smp}_final.vcf.gz.csi",
           fasta = "db/H37Rv.fna",
           bam = "BAM/{smp}.bam",
           bai = "BAM/{smp}.bam.bai"
    output: "VCF/{smp}_final.duphold.vcf"
    log: "logs/duphold_{smp}.log"
    conda: "envs/duphold.yaml"
    message: """--- Estimating depth of flanking regions in {wildcards.smp}.bam ---"""
    shell: """
        duphold -v {input.vcf} -b {input.bam} -f {input.fasta} -o {output} 2> {log}
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
    log: "logs/survivor_merge.log"
    shell: """
        SURVIVOR merge {input} 5000 1 1 0 0 100 {output} 2> {log}
  """

## Annotate VCF files
rule annotate_vcf:
    input: 'VCF/merged_raw.vcf'
    output: "VCF/merged.vcf"
    message: """--- Annotating VCF ---"""
    log: "logs/snpeff.log"
    shell: """
        snpEff ann m_tuberculosis_H37Rv -noLog -noStats -no-downstream -no-upstream -no-utr -download -o vcf {input} > {output} 2> {log}
  """

## Convert multisample VCF to tab-separated format
rule convert2table:
    input: "VCF/merged.vcf"
    output: "RD_table_raw.tsv"
    message: """--- Converting VCF to tab-separated table ---"""
    log: "logs/gatk_vcf2table.log"
    shell: """
        gatk VariantsToTable -V {input} -F CHROM -F POS -F END -F SVLEN -F SVTYPE -F LOF -GF LN -O {output} 2> {log}
  """

## Remove unnecessary info from the table
rule clean_table:
    input: "RD_table_raw.tsv"
    output: "RD_table_cleaned.tsv"
    message: """--- Cleaning table ---"""
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
    message: """--- Annotating putative RDs ---"""
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
    message: """--- Computing deleted regions coverage in {wildcards.smp}.bam ---"""
    shell: """
        parallel mosdepth -n -m -x -b {input.rd} {params.prefix} ::: {input.bam}
  """

## Decompress BED files
rule decompress:
    input: bedgz = "BED/{smp}.regions.bed.gz"
    output: bed = "BED/{smp}.regions.bed"
    threads: 2
    message: """--- Decompressing {wildcards.smp}.regions.bed.gz ---"""
    shell: """
        pigz -p {threads} -dc {input.bedgz} > {output.bed}
  """

## Calculate the proportion of the read depth in RD regions to the total chromosome depth
rule calculate_proportion:
    input: "BED/{smp}.regions.bed"
    output: "BED/{smp}.proportion.bed"
    message: """--- Calculating read depth ratio of RD regions in {wildcards.smp} ---"""
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
    message: """--- Concatenating BED files ---"""
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
    message: """--- Making final tables ---"""
    log: "logs/r_tables.log"
    shell: "Rscript scripts/makeTables.R -T {config[threshold]} -i {input} -x {output.xlsx} -t {output.tab} -b {output.bin} 2> {log}"

onstart:
    date_st = datetime.datetime.now().strftime("%m-%d-%Y")
    time_st = datetime.datetime.now().strftime("%H:%M:%S")
    print("Running RDscan pipeline...")
    print("Date: " + date_st + " " + "Time: " + time_st)

onsuccess:
    dateStr = datetime.datetime.now().strftime("%m-%d-%Y_%Hh%Mm%Ss")
    date_ss = datetime.datetime.now().strftime("%m-%d-%Y")
    time_ss = datetime.datetime.now().strftime("%H:%M:%S")
    t_end = time.time()
    hours, rem = divmod(t_end-t_start, 3600)
    minutes, seconds = divmod(rem, 60)
    resultsDir = "Results" + "_" + dateStr
    os.makedirs(resultsDir, exist_ok=True)
    shutil.move("RD_known.bin.tsv", resultsDir)
    shutil.move("RD_known.tsv", resultsDir)
    shutil.move("RD_known.xlsx", resultsDir)
    shutil.move("RD_table.tsv", resultsDir)
    # shutil.rmtree(".snakemake", ignore_errors=True)
    shutil.rmtree("VCF", ignore_errors=True)
    shutil.rmtree("BED", ignore_errors=True)
    os.remove("RD_table_raw.tsv")
    os.remove("RD_table_cleaned.tsv")
    os.remove("all_concat.bed")
    print("Completed successfully!")
    print("Elapsed time: " + "{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))
    print("Date: " + date_ss + " " + "Time: " + time_ss)

onerror:
    date_er = datetime.datetime.now().strftime("%m-%d-%Y")
    time_er = datetime.datetime.now().strftime("%H:%M:%S")
    print("An error occurred!")
    print("Date: " + date_er + " " + "Time: " + time_er)
