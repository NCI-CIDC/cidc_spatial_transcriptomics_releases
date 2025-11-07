rule azimuth:
    input:
        done=rules.install_azimuth_ref.output.done,
        rds=rules.cell_cycle_scoring.output.rds if config["cell_cycle_scoring"]["cell_cycle_score"]=="TRUE" else rules.dimensionality_reduction.output.rds
    output:
        txt=paths.cell_type_annotation.txt_azimuth,
        rds=paths.cell_type_annotation.rds_azimuth
    benchmark:
        "benchmark/{sample}_azimuth.tab"
    log:
        "log/{sample}_azimuth.log"
    conda:
        "../envs/seurat.yaml"
    params:
        ref=config["cell_type_annotation"]["azimuth_ref"],
        lvl=config["cell_type_annotation"]["azimuth_lvl"]
    threads: config["azimuth_threads"]
    script:
        "../source/r/azimuth.R"

rule tabula_sapiens:
    input:
        ref=rules.retrieve_ts_ref.output.ts_ref,
        rds=TS_INPUT
    output:
        txt=paths.cell_type_annotation.txt_ts,
        rds=paths.cell_type_annotation.rds_ts
    benchmark:
        "benchmark/{sample}_tabula_sapiens.tab"
    log:
        "log/{sample}_tabula_sapiens.log"
    conda:
        "../envs/seurat.yaml"
    params:
        sctransform=config["normalization"]["sctransform"],
        lvl=config["cell_type_annotation"]["ts_lvl"],
        k_anchor=config["cell_type_annotation"]["ts_k_anchor"]
    threads: config["ts_threads"]
    script:
        "../source/r/tabula-sapiens.R"

## Input was originally ref=rules.retrieve_custom_ref.output.custom_ref
rule custom:
    input:
        ref=rules.create_custom_ref.output.rds,
        rds=CUSTOM_INPUT,
        txt=paths.ref_files.custom_txt
    output:
        txt=paths.cell_type_annotation.txt_custom,
        rds=paths.cell_type_annotation.rds_custom
    benchmark:
        "benchmark/{sample}_custom.tab"
    log:
        "log/{sample}_custom.log"
    conda:
        "../envs/seurat.yaml"
    params:
        sctransform=config["normalization"]["sctransform"]
    threads: config["custom_threads"]
    script:
        "../source/r/custom.R"
