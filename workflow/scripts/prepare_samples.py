import sys
import pandas as pd

# Redirect stderr to the Snakemake log file
log_fh = open(snakemake.log[0], "w")
sys.stdout = log_fh
sys.stderr = log_fh

# Read inputs, outputs and parameters from the Snakemake object
counts_path = snakemake.input.counts
metadata_path = snakemake.input.metadata
out_path = snakemake.output[0]
condition_column = snakemake.params.condition_column

# Load raw counts and sample metadata
counts_df = pd.read_csv(counts_path, index_col=0, low_memory=False)
metadata = pd.read_csv(metadata_path, index_col=0, low_memory=False)

# Align sample metadata to match count matrix columns
common_samples = counts_df.columns.intersection(metadata.index)
metadata = metadata.loc[common_samples].copy()

# Ensure condition column exists
if condition_column not in metadata.columns:
    raise ValueError(f"Column {condition_column} not found in metadata")

# Save aligned sample metadata table in samples.tsv
metadata.to_csv(out_path, sep="\t", index=True, index_label="sample_id")