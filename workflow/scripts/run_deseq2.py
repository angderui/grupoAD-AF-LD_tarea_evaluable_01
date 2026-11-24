import sys
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Redirect stdout and stderr to the Snakemake log file
log_fh = open(snakemake.log[0], "w")
sys.stdout = log_fh
sys.stderr = log_fh

# Read inputs, outputs, and parameters from Snakemake
samples_path = snakemake.input.samples
counts_path = snakemake.input.counts
results_path = snakemake.output.results

condition_column = snakemake.params.condition_column
reference = snakemake.params.reference
test = snakemake.params.test
covariates = list(snakemake.params.covariates) if snakemake.params.covariates else []
min_residual_df = snakemake.params.min_residual_df
min_count = snakemake.params.min_count
alpha = snakemake.params.alpha

subset_name = snakemake.wildcards.subset

print(f"[INFO] Starting DESeq2 for subset '{subset_name}'")
print(f"[INFO] Reading samples from: {samples_path}")
print(f"[INFO] Reading counts from: {counts_path}")

samples_df = pd.read_csv(samples_path, sep="\t", index_col=0)
counts_df = pd.read_csv(counts_path, index_col=0, low_memory=False)

# Align sample IDs between count matrix and sample metadata
common_samples = samples_df.index.intersection(counts_df.columns)
if len(common_samples) == 0:
    raise ValueError("No overlapping sample IDs between samples table and count matrix.")

samples_df = samples_df.loc[common_samples].copy()
counts_df = counts_df.loc[:, common_samples].copy()

assert list(counts_df.columns) == list(samples_df.index), "Sample alignment mismatch."

print(f"[INFO] Aligned samples: {len(common_samples)} | Initial genes: {counts_df.shape[0]}")

# Validate condition column and set reference order
if condition_column not in samples_df.columns:
    raise ValueError(f"Condition column '{condition_column}' not found in samples table.")

samples_df[condition_column] = pd.Categorical(
    samples_df[condition_column],
    categories=[reference, test],
    ordered=True,
)
if samples_df[condition_column].isna().any():
    raise ValueError(
        f"Unexpected or missing values in '{condition_column}'. Expected '{reference}' and '{test}'."
    )

# Filter genes using smallest condition group size in the current subset
group_sizes = samples_df[condition_column].value_counts()
min_samples = int(group_sizes.min())
if min_samples < 1:
    raise ValueError(f"Invalid group sizes in subset '{subset_name}': {group_sizes.to_dict()}")

keep_genes = (counts_df >= min_count).sum(axis=1) >= min_samples
counts_filtered_df = counts_df.loc[keep_genes].copy()

print(
    f"[INFO] Group sizes ({condition_column}): {group_sizes.to_dict()} -> "
    f"min_samples={min_samples}, min_count={min_count}"
)
print(f"[INFO] Filtered genes: {counts_filtered_df.shape[0]} / {counts_df.shape[0]}")

# Decide whether covariates can be retained based on variability and min_residual_df
valid_covariates = []
for cov in covariates:
    if cov not in samples_df.columns:
        print(f"[WARN] Covariate '{cov}' not found in samples table; skipping.")
        continue
    n_levels = samples_df[cov].nunique(dropna=True)
    if n_levels < 2:
        print(f"[WARN] Covariate '{cov}' has fewer than 2 levels ({n_levels}); dropping.")
        continue
    valid_covariates.append(cov)

# Approximate number of model parameters: 1 (intercept) + 1 (condition contrast) + sum(n_levels - 1)
n_total_samples = len(samples_df)
extra_params = sum(int(samples_df[cov].nunique(dropna=True)) - 1 for cov in valid_covariates)
approx_params = 2 + extra_params
approx_residual_df = n_total_samples - approx_params

if valid_covariates and approx_residual_df >= min_residual_df:
    active_covariates = valid_covariates
    for cov in active_covariates:
        samples_df[cov] = samples_df[cov].astype("category")
    design_formula = f"~ {condition_column} + {' + '.join(active_covariates)}"
    print(
        f"[INFO] Using design '{design_formula}' "
        f"(approx_residual_df={approx_residual_df} >= min_residual_df={min_residual_df})"
    )
else:
    active_covariates = []
    design_formula = f"~ {condition_column}"
    print(
        f"[WARN] Falling back to design '{design_formula}' for subset '{subset_name}' "
        f"(covariates={valid_covariates}, approx_residual_df={approx_residual_df}, "
        f"min_residual_df={min_residual_df})"
    )

# PyDESeq2 expects samples as rows and genes as columns
counts_mat_filtered = counts_filtered_df.T

dds = DeseqDataSet(
    counts=counts_mat_filtered,
    metadata=samples_df,
    design=design_formula,
    refit_cooks=True,
)
dds.deseq2()

stat_res = DeseqStats(
    dds,
    contrast=[condition_column, test, reference],
    alpha=alpha,
)
stat_res.summary()

res_df = stat_res.results_df.copy()
res_df.insert(0, "gene_id", res_df.index)
res_df["minus_log10_padj"] = -np.log10(res_df["padj"].fillna(1.0) + 1e-300)
res_df = res_df.sort_values(by="padj", na_position="last")

res_df.to_csv(results_path, sep="\t", index=False)
print(f"[INFO] Saved full DESeq2 results ({len(res_df)} genes) to: {results_path}")