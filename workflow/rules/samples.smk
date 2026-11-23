rule prepare_full:
    input:
        counts=config["data"]["counts"],
        metadata=config["data"]["metadata"],
    output:
        "results/full/samples.tsv",
    log:
        "logs/full/prepare_full.log",
    conda:
        "../envs/python.yaml"
    params:
        condition_column=config["design"]["condition_column"],
    script:
        "../scripts/prepare_samples.py"