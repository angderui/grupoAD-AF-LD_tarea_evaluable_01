import sys
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Redirect stderr to the Snakemake log file
sys.stderr = open(snakemake.log[0], "w")

# Read inputs, outputs and parameters
samples_path = snakemake.input.samples
counts_path = snakemake.input.counts
results_path = snakemake.output.results
de_genes_path = snakemake.output.de_genes

condition_column = snakemake.params.condition_column
reference = snakemake.params.reference
test = snakemake.params.test
covariates = list(snakemake.params.covariates)
min_count = snakemake.params.min_count
alpha = snakemake.params.alpha
lfc_threshold = snakemake.params.lfc_threshold

# Load prepared sample metadata and raw count matrix
samples_df = pd.read_csv(samples_path, sep="\t", index_col=0)
counts_df = pd.read_csv(counts_path, index_col=0, low_memory=False)

# Align count matrix columns to match samples_df index
common_samples = samples_df.index.intersection(counts_df.columns)
samples_df = samples_df.loc[common_samples].copy()
counts_df = counts_df.loc[:, common_samples]

assert list(counts_df.columns) == list(samples_df.index), "Error: Samples do not match."

# Ensure categorical types and reference level order
samples_df[condition_column] = pd.Categorical(
    samples_df[condition_column],
    categories=[reference, test],
    ordered=True,
)
for cov in covariates:
    samples_df[cov] = samples_df[cov].astype("category")

# Filter low_count genes using the smallest condition group size as min_samples
group_sizes = samples_df[condition_column].value_counts()
min_samples = int(group_sizes.min())

keep_genes = (counts_df >= min_count).sum(axis=1) >= min_samples
counts_filtered_df = counts_df.loc[keep_genes]

# PyDESeq2 expects samples in rows and genes in columns
counts_mat_filtered = counts_filtered_df.T

# Construct design formula string (e.g., "~ site + plate")
design_formula = f"~ {condition_column}"
if covariates:
    design_formula += f" + {' + '.join(covariates)}"

# Run PyDESeq2 using the modern design API
dds = DeseqDataSet(
    counts=counts_mat_filtered,
    metadata=samples_df,
    design=design_formula,
    refit_cooks=True,
)
dds.deseq2()

# Extract contrast statistics (test vs reference)
stat_res = DeseqStats(
    dds,
    contrast=[condition_column, test, reference],
    alpha=alpha,
)
stat_res.summary()

res_df = stat_res.results_df.copy()
res_df.insert(0, "gene_id", res_df.index)
res_df["minus_log10_padj"] = -np.log10(res_df["padj"] + 1e-300)

# Sort full results by adjusted p-value before saving
res_df = res_df.sort_values(by="padj", na_position="last")
res_df.to_csv(results_path, sep="\t", index=False)

# Filter differentially expressed genes by padj and |log2FoldChange|
de_mask = (
    res_df["padj"].notna()
    & (res_df["padj"] < alpha)
    & (res_df["log2FoldChange"].abs() > lfc_threshold)
)
de_df = res_df.loc[de_mask].copy()
de_df.to_csv(de_genes_path, sep="\t", index=False)