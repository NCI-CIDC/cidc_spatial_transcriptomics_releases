rule merge:
    input:
        rds=expand(rules.spatial_variable_features.output.rds, sample=SAMID)
    output:
        rds=paths.merge.rds
    benchmark:
        "benchmark/merge.tab"
    log:
        "log/merge.log"
    conda:
        "../envs/seurat.yaml"
    threads: config["merge_threads"]
    script:
        "../source/r/merge.R"
