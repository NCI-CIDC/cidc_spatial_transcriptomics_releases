rule normalization:
    input:
        rds=rules.filter.output.rds
    output:
        rds=paths.normalization.rds,
    benchmark:
        "benchmark/{sample}_normalization.tab"
    log:
        "log/{sample}_normalization.log"
    conda:
        "../envs/seurat.yaml"
    params:
        functions=SOURCEDIR+"/r/functions.R",
        sctransform=config["normalization"]["sctransform"],
        imputation=config["normalization"]["imputation"],
        nfeatures=config["normalization"]["nfeatures"],
        visium_type=lambda wildcards: image_path(wildcards.sample, "visium_type"),
        bin_size=config["qc"]["bin_size"]
    threads: config["normalization_threads"]
    script:
        "../source/r/normalization.R"
