rule check_additional_metadata:
    input:
        metadata=METADATA_CSV
    output:
        csv=paths.qc.csv
    benchmark:
        "benchmark/check_additional_metadata.tab"
    log:
        "log/check_additional_metadata.log"
    script:
        "../source/r/check-additional-metadata.R"

rule create_rds:
    input:
        packages=rules.install_packages.output.csv,
        txt=rules.spaceranger.output.txt,
        csv=rules.check_additional_metadata.output.csv
    output:
        rds=paths.qc.rds
    benchmark:
        "benchmark/{sample}_create_rds.tab"
    log:
        "log/{sample}_create_rds.log"
    conda:
        "../envs/seurat.yaml"
    params:
        outs=paths.spaceranger.outs,
        bin_size=config["qc"]["bin_size"],
        workflow=config["workflow"],
        visium_type=lambda wildcards: image_path(wildcards.sample, "visium_type"),
        sample="{sample}",
    script:
        "../source/r/create-rds.R"

rule qc:
    input:
        rds=rules.create_rds.output.rds
    output:
        plot=paths.qc.plot,
        nfeature=paths.qc.nfeature,
        ncount=paths.qc.ncount,
        mt=paths.qc.mt,
        hb=paths.qc.hb,
        rb=paths.qc.rb,
        nspot=paths.qc.nspot,
        rds=paths.qc.rds_qc
    benchmark:
        "benchmark/{sample}_qc.tab"
    log:
        "log/{sample}_qc.log"
    conda:
        "../envs/seurat.yaml"
    params:
        functions=SOURCEDIR+"/r/functions.R",
        visium_type=lambda wildcards: image_path(wildcards.sample, "visium_type"),
        bin_size=config["qc"]["bin_size"],
    script:
        "../source/r/qc.R"
