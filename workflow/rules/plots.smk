rule plot_deseq2_full:
    input:
        results="results/full/deseq2/results.tsv",
    output:
        volcano="results/full/plots/volcano.png",
        ma="results/full/plots/ma.png",
    log:
        "logs/full/plots.log",
    conda:
        "../envs/python.yaml"
    params:
        alpha=config["deseq2"]["alpha"],
        lfc_threshold=config["deseq2"]["lfc_threshold"],
    script:
        "../scripts/plots.py"
