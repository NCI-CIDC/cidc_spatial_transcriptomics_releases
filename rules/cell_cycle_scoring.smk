rule cell_cycle_scoring:
    input:
        rds=rules.dimensionality_reduction.output.rds
    output:
        rds=paths.cell_cycle_scoring.rds
    benchmark:
        "benchmark/{sample}_cell_cycle_scoring.tab"
    log:
        "log/{sample}_cell_cycle_scoring.log"
    conda:
        "../envs/seurat.yaml"
    threads: config["cell_cycle_scoring_threads"]
    script:
        "../source/r/cell-cycle-scoring.R"
