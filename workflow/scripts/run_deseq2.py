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