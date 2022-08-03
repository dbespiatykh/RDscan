import os
import pandas as pd

df = pd.read_csv(
    snakemake.input[0],
    compression="gzip",
    sep="\t",
    names=["chrom", "start", "end", "RD", "depth"],
    dtype={"chrom": str, "start": str, "end": str, "RD": str, "depth": int},
)

df.loc[:, "depth"] = df.loc[:, "depth"].div(df.iloc[0]["depth"])
df["Sample"] = os.path.basename(snakemake.input[0]).split(".regions.bed.gz")[0]
df = df[df.RD != "H37Rv"]

df.to_csv(snakemake.output[0], sep="\t", index=False)
