rule annotation_postprocessing:
    input:
        rds=ANNOT_POSTPROCESS_INPUT
    output:
        rds=paths.annotation_postprocessing.rds,
        txt=paths.annotation_postprocessing.txt,
        png=paths.annotation_postprocessing.png if config["annotation_postprocessing"]["filter"] == "TRUE" else [],
    benchmark:
        "benchmark/{sample}_annotation_postprocessing.tab"
    log:
        "log/{sample}_annotation_postprocessing.log"
    conda:
        "../envs/seurat.yaml"
    params:
        subset=config["annotation_postprocessing"]["subset"],
        ref=config["annotation_postprocessing"]["ref"],
        azimuth_ref=config["cell_type_annotation"]["azimuth_ref"],
        ts_ref=config["cell_type_annotation"]["ts_ref"],
        lvl=config["annotation_postprocessing"]["lvl"],
        type=config["annotation_postprocessing"]["type"],
        filter=config["annotation_postprocessing"]["filter"],
        filter_ref=config["annotation_postprocessing"]["filter_ref"],
        filter_lvl=config["annotation_postprocessing"]["filter_lvl"],
        filter_score=config["annotation_postprocessing"]["filter_score"],
        functions=SOURCEDIR+"/r/functions.R"
    script:
        "../source/r/annotation-postprocessing.R"
