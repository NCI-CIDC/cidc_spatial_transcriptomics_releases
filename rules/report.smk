rule render_report:
    input:
        rds=rules.spatial_variable_features.output.rds,
        metadata=rules.check_additional_metadata.output.csv,
        qc_plot_pre=rules.qc.output.plot,
        qc_nspot_pre=rules.qc.output.nspot,
        ncount_pre=rules.qc.output.ncount,
        mt_pre=rules.qc.output.mt,
        qc_plot_post=rules.filter.output.plot,
        ncount_post=rules.filter.output.ncount,
        mt_post=rules.filter.output.mt,
        azimuth_txt=rules.azimuth.output.txt if config["cell_type_annotation"]["azimuth"] == "TRUE" else [],
        ts_txt=rules.tabula_sapiens.output.txt if config["cell_type_annotation"]["ts"] == "TRUE" else [],
        custom_txt=rules.custom.output.txt if config["cell_type_annotation"]["custom"] == "TRUE" else [],
        cluster_umap=rules.clustering.output.umap,
        cluster_individual=rules.clustering.output.individual,
        svg_csv=rules.spatial_variable_features.output.csv
    output:
        heatmap=paths.report.heatmap,
        azimuth_a=paths.report.azimuth_a if config["cell_type_annotation"]["azimuth"] == "TRUE" else [],
        azimuth_aa=paths.report.azimuth_aa if config["cell_type_annotation"]["azimuth"] == "TRUE" else [],
        ts_a=paths.report.ts_a if config["cell_type_annotation"]["ts"] == "TRUE" else [],
        ts_aa=paths.report.ts_aa if config["cell_type_annotation"]["ts"] == "TRUE" else [],
        custom_a=paths.report.custom_a if config["cell_type_annotation"]["custom"] == "TRUE" else [],
        custom_aa=paths.report.custom_aa if config["cell_type_annotation"]["custom"] == "TRUE" else [],
        r=paths.report.r,
        html=paths.report.html
    benchmark:
        "benchmark/{sample}_render_report.tab"
    log:
        "log/{sample}_render_report.log"
    conda:
        "../envs/seurat.yaml"
    params:
       functions=SOURCEDIR+"/r/functions.R",
       dir=PREDIR+"/report",
       template=SOURCEDIR+"/r/report-template.R",
       sample="{sample}",
       sctransform=config["normalization"]["sctransform"],
       hla_ig=config["filter"]["hla_ig"],
       mt=config["filter"]["mt"],
       rb=config["filter"]["rb"],
       hb=config["filter"]["hb"],
       ncount_min=config["filter"]["ncount_min"],
       nfeature_min=config["filter"]["nfeature_min"],
       visium_type=lambda wildcards: image_path(wildcards.sample, "visium_type"),
       bin_size=config["qc"]["bin_size"],
       azimuth=config["cell_type_annotation"]["azimuth"],
       azimuth_ref=config["cell_type_annotation"]["azimuth_ref"],
       azimuth_lvl=config["cell_type_annotation"]["azimuth_lvl"],
       ts=config["cell_type_annotation"]["ts"],
       ts_lvl=config["cell_type_annotation"]["ts_lvl"],
       custom=config["cell_type_annotation"]["custom"],
       azimuth_b=paths.report.azimuth_b,
       azimuth_bb=paths.report.azimuth_bb,
       ts_b=paths.report.ts_b,
       ts_bb=paths.report.ts_bb,
       custom_b=paths.report.custom_b,
       custom_bb=paths.report.custom_bb,
       svg_plot=paths.spatial_variable_features.png
    threads: config["report_threads"]
    script:
        "../source/r/render-report.R"
