from snakemake.utils import min_version

min_version("9.0")


SIZES = config["subsampling"]["sizes"]
SEEDS = config["subsampling"]["seeds"]

# One analysis per (size, seed) pair, plus the full cohort
SUBSETS = [f"n{n}_s{seed}" for n in SIZES for seed in SEEDS]
ALL_SUBSETS = ["full"] + SUBSETS

LIBRARIES = list(config["enrichment"]["libraries"])
ENRICH_SUBSETS = config["enrichment"]["subsets"]
DIRECTIONS = ["up", "down"]


# "full" or "n<size>_s<seed>"; nothing else is a valid subset
wildcard_constraints:
    subset="full|n[0-9]+_s[0-9]+",
    n="[0-9]+",    
    seed="[0-9]+",
    library="|".join(LIBRARIES),
    direction="up|down",