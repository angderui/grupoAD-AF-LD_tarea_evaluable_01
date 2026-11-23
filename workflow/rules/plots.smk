rule plot_volcano_full:
    input:
        results="results/full/deseq2/results.tsv",
    output:
        volcano="results/full/plots/volcano.png",
    log:
        "logs/full/plot_volcano.log",
    conda:
        "../envs/python.yaml"
    params:
        alpha=config["deseq2"]["alpha"],
        lfc_threshold=config["deseq2"]["lfc_threshold"],
    script:
        "../scripts/plot_volcano.py"

rule plot_ma_full:
    input:
        results="results/full/deseq2/results.tsv",
    output:
        ma="results/full/plots/ma.png",
    log:
        "logs/full/plot_ma.log",
    conda:
        "../envs/python.yaml"
    params:
        alpha=config["deseq2"]["alpha"],
        lfc_threshold=config["deseq2"]["lfc_threshold"],
    script:
        "../scripts/plot_ma.py"    
