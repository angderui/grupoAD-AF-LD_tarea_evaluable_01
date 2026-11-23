import sys
import matplotlib

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
df = pd.read_csv(results_path, sep="\t")

# Ensure required columnns exists

# Volcano plot
