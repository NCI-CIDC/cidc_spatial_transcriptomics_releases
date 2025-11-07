rule clustering:
    input:
        rds=CLUSTERING_INPUT
    output:
        umap=paths.clustering.umap,
        individual=paths.clustering.individual,
        rds=paths.clustering.rds
    benchmark:
        "benchmark/{sample}_clustering.tab"
    log:
        "log/{sample}_clustering.log"
    conda:
        "../envs/seurat.yaml"
    params:
        samid=SAMID,
        visium_type=lambda wildcards: image_path(wildcards.sample, "visium_type")
    threads: config["clustering_threads"]
    script:
        "../source/r/clustering.R"
