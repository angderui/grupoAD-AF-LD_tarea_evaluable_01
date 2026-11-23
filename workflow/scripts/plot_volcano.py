import sys
import matplotlib

# Use non-interactive backend so plots can be saved without a display
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Redirect stdout and stderr to the Snakemake log file
log_fh = open(snakemake.log[0], "w")
sys.stdout = log_fh
sys.stderr = log_fh

# Read inputs, outputs and parameters
results_path = snakemake.input.results
volcano_out = snakemake.output.volcano

alpha = snakemake.params.alpha
lfc_threshold = snakemake.params.lfc_threshold

# Load results.tsv
print(f"Reading DESeq2 results from: {results_path}")
df = pd.read_csv(results_path, sep="\t")

# Ensure gene_id column is present
if "gene_id" not in res_df.columns:
    res_df["gene_id"] = res_df.index.astype(str)

# Ensure -log10(padj) column is present
if "minus_log10_padj" not in res_df.columns:
    res_df["minus_log10_padj"] = -np.log10(res_df["padj"].fillna(1.0) + 1e-300)

# Prepare valid rows for Volcano plot
res_plot = res_df.dropna(subset=["padj", "log2FoldChange", "minus_log10_padj"]).copy()

# Define masks to separate significant Up/Down genes from non-significant genes
up_mask = (res_plot["padj"] < alpha) & (res_plot["log2FoldChange"] > lfc_threshold)
down_mask = (res_plot["padj"] < alpha) & (res_plot["log2FoldChange"] < -lfc_threshold)
ns_mask = ~(up_mask | down_mask)

# Select top 10 Up and top 10 Down genes by adjusted p-value
top10_up = (
    res_plot.loc[up_mask]
    .sort_values(by=["padj", "log2FoldChange"], ascending=[True, False])
    .head(10)
)
top10_down = (
    res_plot.loc[down_mask]
    .sort_values(by=["padj", "log2FoldChange"], ascending=[True, True])
    .head(10)
)
top10_volcano = pd.concat([top10_up, top10_down])

# Create Volcano plot figure
plt.figure(figsize=(8, 6))

# Plot each group separately with its label for the legend
plt.scatter(
    res_plot.loc[ns_mask, "log2FoldChange"],
    res_plot.loc[ns_mask, "minus_log10_padj"],
    c="darkgray",
    s=10,
    alpha=0.5,
    label="No significativo",
)
plt.scatter(
    res_plot.loc[up_mask, "log2FoldChange"],
    res_plot.loc[up_mask, "minus_log10_padj"],
    c="firebrick",
    s=10,
    alpha=0.7,
    label="Up",
)
plt.scatter(
    res_plot.loc[down_mask, "log2FoldChange"],
    res_plot.loc[down_mask, "minus_log10_padj"],
    c="steelblue",
    s=10,
    alpha=0.7,
    label="Down",
)

# Draw threshold lines
plt.axvline(lfc_threshold, color="black", linestyle="--", alpha=0.6)
plt.axvline(-lfc_threshold, color="black", linestyle="--", alpha=0.6)
plt.axhline(-np.log10(alpha), color="black", linestyle="--", alpha=0.6)

# Add text labels for top Up and Down genes
for _, row in top10_volcano.iterrows():
    plt.text(
        x=row["log2FoldChange"],
        y=row["minus_log10_padj"],
        s=row["gene_id"],
        fontsize=8,
        weight="bold",
        ha="right",
        va="bottom",
    )

# Set axis labels, title, and legend
plt.xlabel("log2FoldChange")
plt.ylabel("-log10(padj)")
plt.title("Volcano Plot (Tumor vs Normal)")
plt.legend(loc="upper left", frameon=True)

# Save figure and close plot
plt.tight_layout()
plt.savefig(volcano_out, dpi=300)
plt.close()

print(f"Volcano plot saved successfully to: {volcano_out}")