import sys
import matplotlib

# Use non-interactive backend so plots can be saved without a display
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import gseapy as gp

# Redirect stdout and stderr to the Snakemake log file
log_fh = open(snakemake.log[0], "w")
sys.stdout = log_fh
sys.stderr = log_fh

# Read inputs, outputs and parameters
de_path = snakemake.input.de_genes
results_path = snakemake.input.results
table_out = snakemake.output.table
plot_out = snakemake.output.plot

subset = snakemake.wildcards.subset
library_key = snakemake.wildcards.library
direction = snakemake.wildcards.direction
library_name = snakemake.params.library_name
max_term = int(snakemake.params.max_term)

print(f"[INFO] Subset={subset} | library={library_key} ({library_name}) | direction={direction}")
print(f"[INFO] Reading DE genes from: {de_path}")
print(f"[INFO] Reading background results from: {results_path}")

# Load differentially expressed genes and full DESeq2 results (used as background universe)
de_df = pd.read_csv(de_path, sep="\t")
res_all_df = pd.read_csv(results_path, sep="\t")

# Validate required columns before filtering or querying Enrichr
if "gene_id" not in de_df.columns or "log2FoldChange" not in de_df.columns:
    raise ValueError(f"Missing 'gene_id' or 'log2FoldChange' in {de_path}")
if "gene_id" not in res_all_df.columns:
    raise ValueError(f"Missing 'gene_id' in {results_path}")

# Extract unique tested genes to define the enrichment background
background_genes = res_all_df["gene_id"].dropna().astype(str).unique().tolist()

# Select up- or down-regulated DE genes according to the requested wildcard direction
if direction == "up":
    gene_list = (
        de_df.loc[de_df["log2FoldChange"] > 0, "gene_id"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )
elif direction == "down":
    gene_list = (
        de_df.loc[de_df["log2FoldChange"] < 0, "gene_id"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )
else:
    raise ValueError(f"Unsupported direction: {direction}")

print(f"[INFO] Target genes ({direction}): {len(gene_list)} | Background genes: {len(background_genes)}")

# Handle empty gene lists gracefully by writing an empty table and placeholder figure
if len(gene_list) == 0:
    print(f"[WARN] No {direction} DE genes in '{subset}'. Writing empty table and placeholder plot.")
    pd.DataFrame().to_csv(table_out, sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5, f"No {direction} DE genes ({subset})", ha="center", va="center")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(plot_out, dpi=300)
    plt.close(fig)
    sys.exit(0)

# Query Enrichr via gseapy using the tested gene universe as background
enr = gp.enrichr(
    gene_list=gene_list,
    gene_sets=library_name,
    organism="human",
    background=background_genes,
    outdir=None,
    no_plot=True,
)

# Save full enrichment results table to TSV
df_res = enr.res2d.copy() if enr.res2d is not None else pd.DataFrame()
df_res.to_csv(table_out, sep="\t", index=False)
print(f"[INFO] Saved enrichment table ({len(df_res)} terms) to: {table_out}")

# Generate a dotplot for top enriched terms, or a placeholder figure if no terms are returned
if not df_res.empty:
    gp.dotplot(
        df_res,
        title=f"{subset} - {library_key} ({direction})",
        cmap="viridis",
        max_term=max_term,
        ofname=plot_out,
    )
    plt.close("all")
    print(f"[INFO] Saved dotplot to: {plot_out}")
else:
    print(f"[WARN] Enrichr returned empty results for {subset} {library_key} ({direction}).")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5, f"No enriched terms ({subset}, {library_key}, {direction})", ha="center", va="center")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(plot_out, dpi=300)
    plt.close(fig)