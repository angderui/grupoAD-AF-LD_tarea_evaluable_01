rule run_deseq2:
    input:
        samples="results/{subset}/samples.tsv",
        counts=config["data"]["counts"],
    output:
        results="results/{subset}/deseq2/results.tsv",
    log:
        "logs/{subset}/run_deseq2.log",
    conda:
        "../envs/pydeseq2.yaml"
    params:
        condition_column=config["design"]["condition_column"],
        reference=config["design"]["reference"],
        test=config["design"]["test"],
        covariates=config["design"]["covariates"],
        min_residual_df=config["design"]["min_residual_df"],
        min_count=config["filter"]["min_count"],
        alpha=config["deseq2"]["alpha"],
    script:
        "../scripts/run_deseq2.py"


rule filter_de_genes:
    input:
        results="results/{subset}/deseq2/results.tsv",
    output:
        de_genes="results/{subset}/deseq2/de_genes.tsv",
    log:
        "logs/{subset}/filter_de_genes.log",
    conda:
        "../envs/python.yaml"
    params:
        alpha=config["deseq2"]["alpha"],
        lfc_threshold=config["deseq2"]["lfc_threshold"],
    script:
        "../scripts/filter_de_genes.py"