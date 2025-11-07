rule dimensionality_reduction:
    input:
        rds=rules.normalization.output.rds
    output:
        pca=paths.dimensionality_reduction.pca,
        vizdim=paths.dimensionality_reduction.vizdim,
        heatmap=paths.dimensionality_reduction.heatmap,
        elbow=paths.dimensionality_reduction.elbow,
        rds=paths.dimensionality_reduction.rds
    benchmark:
        "benchmark/{sample}_dimensionality_reduction.tab"
    log:
        "log/{sample}_dimensionality_reduction.log"
    conda:
        "../envs/seurat.yaml"
    threads: config["dimensionality_reduction_threads"]
    script:
        "../source/r/dimensionality-reduction.R"
