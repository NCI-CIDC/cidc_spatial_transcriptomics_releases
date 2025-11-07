## Set up directory structure based on dirs supplied in config
## Ignores non-zero exit status returned when any directories already exist
rule directory_setup:
    output:
        dirs="progress/dirs.done"
    params:
        subdirs=SUBDIRS
    threads:1
    shell:
        '''
          mkdir {params.subdirs} -p 2> /dev/null
          touch {output}
        '''

## Download reference genome files
rule retrieve_reference_files:
    input:
        dirs=rules.directory_setup.output.dirs
    output:
        fa=paths.ref_files.fa,
        fai=paths.ref_files.fai,
        gtf=paths.ref_files.gtf,
        refdata=directory(paths.ref_files.refdata)
    benchmark:
        "benchmark/retrieve_reference_files.tab"
    log:
        "log/retrieve_reference_files.log"
    params:
        ref_genome=REF_GENOME,
        cloud=CLOUD
    priority: 1000
    shell:
        '''
          echo "{params.cloud} cp -r {params.ref_genome} ref_files" | tee {log}
          {params.cloud} cp -r {params.ref_genome} ref_files 2>> {log}
        '''

## Download human transcriptome probe sets (version 1 and 2)
rule retrieve_probe_sets:
    input:
        dirs=rules.directory_setup.output.dirs
    output:
        probe_v1=paths.ref_files.probe_v1,
        probe_v2=paths.ref_files.probe_v2
    benchmark:
        "benchmark/retrieve_probe_sets.tab"
    log:
        "log/retrieve_probe_sets.log"
    params:
        probe_v1=PROBE_V1,
        probe_v2=PROBE_V2,
        cloud=CLOUD
    priority: 1000
    shell:
        '''
          echo "{params.cloud} cp {params.probe_v1} {params.probe_v2} ref_files" | tee {log}
          {params.cloud} cp {params.probe_v1} {params.probe_v2} ref_files 2>> {log}
        '''

## Retrieve Space Ranger tar file
rule retrieve_spaceranger:
    input:
        refdata=rules.retrieve_reference_files.output.refdata
    output:
        sr_tar=paths.ref_files.sr_tar,
    benchmark:
        "benchmark/retrieve_spaceranger.tab"
    log:
        "log/retrieve_spaceranger.log"
    params:
        sr_tar=SR_TAR,
        cloud=CLOUD
    shell:
        '''
          echo "{params.cloud} cp {params.sr_tar} ref_files" | tee {log}
          {params.cloud} cp {params.sr_tar} ref_files 2>> {log}
        '''
## Untar Space Ranger
rule set_up_spaceranger:
    input:
        sr_tar=rules.retrieve_spaceranger.output.sr_tar
    output:
        done=PREDIR+"/progress/spaceranger_env.done",
        spaceranger=directory(PREDIR+"/spaceranger-3.1.3")
    benchmark:
        "benchmark/set_up_spaceranger.tab"
    log:
        "log/set_up_spaceranger.log"
    params:
        sourcedir=SOURCEDIR
    shell:
        '''
            echo "tar -zxvf {input.sr_tar} && touch {output.done}" | tee {log}
            tar -zxvf {input.sr_tar} && touch {output.done} 2>> {log}
        '''

## Install SeuratWrappers, DoubletFinder, and SeuratData (version 0.2.2.9001)
rule install_packages:
    output:
        csv=PREDIR+"/R_package_versions.csv"
    benchmark:
        "benchmark/install_packages.tab"
    log:
        "log/install_packages.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../source/r/install-packages.R"

## Install the specified reference for use in rule azimuth for cell type annotation
rule install_azimuth_ref:
    input:
        packages=rules.install_packages.output.csv
    output:
        done=PREDIR+"/progress/install_azimuth_ref.done"
    benchmark:
        "benchmark/install_azimuth_ref.tab"
    log:
        "log/install_azimuth_ref.log"
    conda:
        "../envs/seurat.yaml"
    params:
        ref=config["cell_type_annotation"]["azimuth_ref"]
    script:
        "../source/r/install-azimuth-ref.R"

### Retrieve the specified reference for use in rule tabula_sapiens for cell type annotation
rule retrieve_ts_ref:
    output:
        ts_ref=TS_PATH
    benchmark:
        "benchmark/retrieve_ts_ref.tab"
    log:
        "log/retrieve_ts_ref.log"
    params:
        uri=config["cell_type_annotation"]["ts_ref"],
        cloud=CLOUD
    shell:
        '''
          echo "{params.cloud} cp {params.uri} ref_files" | tee {log}
          {params.cloud} cp {params.uri} ref_files 2>> {log}
        '''

### THIS RULE MIGHT BE DELETED!
### Retrieve the specified reference for use in rule custom for cell type annotation
#rule retrieve_custom_ref:
#    output:
#        custom_ref=CUSTOM_PATH
#    benchmark:
#        "benchmark/retrieve_custom_ref.tab"
#    log:
#        "log/retrieve_custom_ref.log"
#    params:
#        uri=config["cell_type_annotation"]["custom_ref"],
#        cloud=CLOUD
#    shell:
#        '''
#          echo "{params.cloud} cp {params.uri} ref_files" | tee {log}
#          {params.cloud} cp {params.uri} ref_files 2>> {log}
#        '''

## Retrieve the counts CSV and metadata CSV to generate the custom reference for cell type annotation
rule retrieve_custom_files:
    output:
        counts=CUSTOM_COUNTS,
        metadata=CUSTOM_METADATA
    benchmark:
        "benchmark/retrieve_custom_files.tab"
    log:
        "log/retrieve_custom_files.log"
    params:
        counts=config["cell_type_annotation"]["custom_counts"],
        metadata=config["cell_type_annotation"]["custom_metadata"],
        cloud=CLOUD
    shell:
        '''
          echo "{params.cloud} cp {params.counts} {params.metadata} ref_files" | tee {log}
          {params.cloud} cp {params.counts} {params.metadata} ref_files 2>> {log}
        '''
## Create custom reference
rule create_custom_ref:
    input:
        counts=rules.retrieve_custom_files.output.counts,
        metadata=rules.retrieve_custom_files.output.metadata,
        packages=rules.install_packages.output.csv
    output:
        txt=paths.ref_files.custom_txt,
        rds=paths.ref_files.custom_ref
    benchmark:
        "benchmark/create_custom_ref.tab"
    log:
        "log/create_custom_ref.log"
    conda:
        "../envs/seurat.yaml"
    threads: config["custom_ref_threads"]
    script:
        "../source/r/create-custom-ref.R"
