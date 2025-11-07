rule spatial_variable_features:
    input:
        rds=rules.clustering.output.rds
    output:
        csv=paths.spatial_variable_features.csv,
        rds=paths.spatial_variable_features.rds
    benchmark:
        "benchmark/{sample}_spatial_variable_features.tab"
    log:
        "log/{sample}_spatial_variable_features.log"
    conda:
        "../envs/seurat.yaml"
    params:
        sctransform=config["normalization"]["sctransform"],
        png=paths.spatial_variable_features.png,
        visium_type=lambda wildcards: image_path(wildcards.sample, "visium_type"),
        sample="{sample}",
        path=PREDIR+"/spatial_variable_features/"
    threads: config["spatial_variable_features_threads"]
    script:
        "../source/r/spatial-variable-features.R"
