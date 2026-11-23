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
ma_out = snakemake.output.ma

alpha = snakemake.params.alpha
lfc_threshold = snakemake.params.lfc_threshold

# Load results.tsv
print(f"Reading DESeq2 results from: {results_path}")
res_df = pd.read_csv(results_path, sep="\t")

# Prepare valid rows for MA plot
res_plot = res_df.dropna(subset=["padj", "log2FoldChange", "baseMean"]).copy()

# Assign point colors by significance and fold change direction
colors = np.where(
    (res_plot["padj"] < alpha) & (res_plot["log2FoldChange"] > lfc_threshold),
    "firebrick",
    np.where(
        (res_plot["padj"] < alpha) & (res_plot["log2FoldChange"] < -lfc_threshold),
        "steelblue",
        "darkgray",
    ),
)

# Create MA plot figure
plt.figure(figsize=(8, 6))
plt.scatter(res_plot["baseMean"], res_plot["log2FoldChange"], c=colors, s=10, alpha=0.7)
plt.axhline(0, color="black", linewidth=1)
plt.xscale("log")
plt.xlabel("mean of normalized counts")
plt.ylabel("log2 fold change")
plt.title("DESeq2 MA Plot")

# Save figure and close plot
plt.tight_layout()
plt.savefig(ma_out, dpi=300)
plt.close()

print(f"MA plot saved successfully to: {ma_out}")