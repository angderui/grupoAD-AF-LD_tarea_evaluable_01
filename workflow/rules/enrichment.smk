rule run_enrichment:
    input:
        de_genes="results/{subset}/deseq2/de_genes.tsv",
        results="results/{subset}/deseq2/results.tsv",
    output:
        table="results/{subset}/enrichment/{library}_{direction}.tsv",
        plot="results/{subset}/enrichment/{library}_{direction}.png",
    log:
        "logs/{subset}/enrichment_{library}_{direction}.log",
    conda:
        "../envs/gseapy.yaml"
    retries: 3
    params:
        library_name=lambda wildcards: config["enrichment"]["libraries"][wildcards.library],
        max_term=config["enrichment"]["top_terms"],
    script:
        "../scripts/run_enrichment.py"