import pandas as pd


def load_data(bed):
    f = pd.read_csv(bed, sep="\t")
    return f

df = pd.concat(load_data(bed) for bed in snakemake.input)

df = df.drop(["chrom", "start", "end"], axis=1)

df.to_csv(snakemake.output[0], sep="\t", index=False)
