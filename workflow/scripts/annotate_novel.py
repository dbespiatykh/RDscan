import numpy as np
import pandas as pd

df1 = pd.read_csv(snakemake.input[0], sep="\t")
df2 = pd.read_csv(snakemake.input[1], sep="\t", names=["CHROM", "START", "END", "RD"])

df2 = df2.iloc[1:-2, :].reset_index(drop=True)
df2["SIZE"] = df2["END"] - df2["START"]

df1 = df1.rename(columns={"POS": "START", "SVLEN": "SIZE", "SVTYPE": "TYPE"})

df1["SIZE"] = df1["SIZE"].abs()
df1.iloc[:, 5:] = df1.iloc[:, 5:].replace({0: np.nan})
df1 = df1[df1.CHROM != "RvD1"]
df1 = df1[df1.CHROM != "TbD1"]
df1.columns = df1.columns.str.replace(r".LN", "", regex=False)

s1 = df1["START"].values
e1 = df1["END"].values
l1 = df1["SIZE"].values[:, None]
s2 = df2["START"].values
e2 = df2["END"].values
l2 = df2["SIZE"].values

r = (
    pd.DataFrame(
        ((e2 - s1[:, None]) > 0)
        & ((s2 - e1[:, None]) < 0)
        & (l2 / l1 >= 0.8)
        & (l2 / l1 <= 1.4)
    )
    .dot(df2.RD + "|")
    .str[:-1]
)

df1["RD"] = r.values
rd = df1["RD"]
df1.drop(labels=["RD"], axis=1, inplace=True)
df1.insert(5, "RD", rd)

df1.to_csv(snakemake.output[0], sep="\t", index=False)
