import sys
import pandas as pd

# Redirect stdout and stderr to the Snakemake log file
log_fh = open(snakemake.log[0], "w")
sys.stdout = log_fh
sys.stderr = log_fh

results_path = snakemake.input.results
de_genes_path = snakemake.output.de_genes
alpha = snakemake.params.alpha
lfc_threshold = snakemake.params.lfc_threshold

subset_name = snakemake.wildcards.subset

print(f"[INFO] Filtering DE genes for subset '{subset_name}' from: {results_path}")
res_df = pd.read_csv(results_path, sep="\t")

required_cols = {"padj", "log2FoldChange"}
missing_cols = required_cols - set(res_df.columns)
if missing_cols:
    raise ValueError(f"Missing required columns in {results_path}: {sorted(missing_cols)}")

res_df = res_df.sort_values(by="padj", na_position="last")

de_mask = (
    res_df["padj"].notna()
    & (res_df["padj"] < alpha)
    & (res_df["log2FoldChange"].abs() > lfc_threshold)
)
de_df = res_df.loc[de_mask].copy()
de_df.to_csv(de_genes_path, sep="\t", index=False)

n_up = int((de_df["log2FoldChange"] > 0).sum()) if len(de_df) > 0 else 0
n_down = int((de_df["log2FoldChange"] < 0).sum()) if len(de_df) > 0 else 0

print(
    f"[INFO] Saved {len(de_df)} DE genes (Up={n_up}, Down={n_down}) "
    f"with padj < {alpha} and |log2FoldChange| > {lfc_threshold} to: {de_genes_path}"
)