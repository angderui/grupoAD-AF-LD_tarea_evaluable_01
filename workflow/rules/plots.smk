rule plot_volcano:
    input:
        results="results/{subset}/deseq2/results.tsv",
    output:
        volcano="results/{subset}}/plots/volcano.png",
    log:
        "logs/{subset}/plot_volcano.log",
    conda:
        "../envs/python.yaml"
    params:
        alpha=config["deseq2"]["alpha"],
        lfc_threshold=config["deseq2"]["lfc_threshold"],
    script:
        "../scripts/plot_volcano.py"

rule plot_ma:
    input:
        results="results/{subset}/deseq2/results.tsv",
    output:
        ma="results/{subset}/plots/ma.png",
    log:
        "logs/{subset}/plot_ma.log",
    conda:
        "../envs/python.yaml"
    params:
        alpha=config["deseq2"]["alpha"],
        lfc_threshold=config["deseq2"]["lfc_threshold"],
    script:
        "../scripts/plot_ma.py"    
