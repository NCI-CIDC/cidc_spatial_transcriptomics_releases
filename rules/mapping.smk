## Perform alignment and quantification of gene expression using `spaceranger count`
rule spaceranger:
    input:
        spaceranger=rules.set_up_spaceranger.output.spaceranger,
        done=rules.set_up_spaceranger.output.done,
        transcriptome_fa=rules.retrieve_reference_files.output.fa,
        transcriptome_fai=rules.retrieve_reference_files.output.fai,
        transcriptome_dir=rules.retrieve_reference_files.output.refdata,
        fastqs_dir=rules.retrieve_fastqs_dir.output.fastqs_dir,
        cytaimage=lambda wildcards: rules.retrieve_cytaimage.output.done if image_path(wildcards.sample, "cytaimage")[0] != "" else [],
        probe_set=lambda wildcards: rules.retrieve_probe_sets.output.probe_v1 if image_path(wildcards.sample, "probe_set_version") == 1 else rules.retrieve_probe_sets.output.probe_v2,
        loupe_alignment=rules.retrieve_loupe_alignment.output.loupe_json,
        image=lambda wildcards: rules.retrieve_image.output.done if image_path(wildcards.sample, "image")[0] != "" else [],
        darkimage=lambda wildcards: rules.retrieve_darkimage.output.done if image_path(wildcards.sample, "darkimage")[0] != "" else [],
        colorizedimage=lambda wildcards: rules.retrieve_colorizedimage.output.done if image_path(wildcards.sample, "colorizedimage")[0] != "" else [],
        sample_param=rules.generate_sample_param.output.sample_param
    output:
        txt=paths.spaceranger.txt
    benchmark:
        "benchmark/{sample}_spaceranger.tab"
    log:
        "log/{sample}_spaceranger.log"
    params:
        run_id="{sample}",
        output_dir=PREDIR+"/spaceranger/{sample}",
        spaceranger_path=PREDIR+"/spaceranger-3.1.3/bin/spaceranger",
        cytaimage=lambda wildcards: image_path(wildcards.sample, "cytaimage")[1] if image_path(wildcards.sample, "cytaimage")[0] != "" else "NA",
        slide=lambda wildcards: image_path(wildcards.sample, "slide"),
        area=lambda wildcards: image_path(wildcards.sample, "area"),
        dapi_index=lambda wildcards: image_path(wildcards.sample, "dapi_index") if image_path(wildcards.sample, "dapi_index") != "" else "NA",
        image=lambda wildcards: image_path(wildcards.sample, "image")[1] if image_path(wildcards.sample, "image")[0] != "" else "NA",
        darkimage=lambda wildcards: image_path(wildcards.sample, "darkimage")[1] if image_path(wildcards.sample, "darkimage")[0] != "" else "NA",
        colorizedimage=lambda wildcards: image_path(wildcards.sample, "colorizedimage")[1] if image_path(wildcards.sample, "colorizedimage")[0] != "" else "NA",
        localmem=config["spaceranger_localmem"],
        localvmem=config["spaceranger_localvmem"]
    priority: 4
    threads: config["spaceranger_threads"]
    script:
        "../source/python/spaceranger.py"

### Run FASTQC
#rule fastqc:
#    input:
#        rules.run_star.output.bam
#    output:
#        targz=paths.fastqc.targz
#    benchmark:
#        'benchmark/{sample}_fastqc.tab'
#    log:
#        'log/{sample}_fastqc.log'
#    conda:
#        SOURCEDIR+"/../envs/fastqc.yaml"
#    params:
#        sample='{sample}',
#        fq_base='fastqc/{sample}.Aligned.sortedByCoord.out_fastqc',
#        fq_zip='fastqc/{sample}.Aligned.sortedByCoord.out_fastqc.zip',
#        fq_html='fastqc/{sample}.Aligned.sortedByCoord.out_fastqc.html'
#    priority: 1
#    threads: 1
#    shell:
#        '''
#          echo "fastqc {input} -q -o fastqc" > {log}
#          fastqc {input} -q -o fastqc 2>> {log}
#
#          ## unzip, remove zipped results, HTML duplicate, and tarball results
#          unzip -qq {params.fq_zip} -d {params.fq_base} && tar -zcf {output.targz} {params.fq_base} && rm -r {params.fq_zip} {params.fq_html} {params.fq_base}
#
#          ## export rule env details
#          conda env export --no-builds > info/fastqc.info
#        '''
#
### Run RSEQC bam_stat.py
#rule bam_qc:
#    input:
#        bam=rules.run_star.output.bam,
#        idx=rules.index_bam.output.idx
#    output:
#        paths.rseqc.bamqc_txt
#    benchmark:
#        'benchmark/{sample}_bam_qc.tab'
#    log:
#        'log/{sample}_bam_qc.log'
#    conda:
#        SOURCEDIR+"/../envs/rseqc.yaml"
#    params:
#        sample='{sample}'
#    priority: 1
#    threads: 1
#    shell:
#        '''
#          echo "bam_stat.py -i {input.bam} > {output}" | tee {log}
#          bam_stat.py -i {input.bam} > {output} 2>> {log}
#        '''
#
### Run RSEQC read_gc.py
#rule bam_gc:
#    input:
#        bam=rules.run_star.output.bam,
#        idx=rules.index_bam.output.idx
#    output:
#        r=paths.rseqc.bamgc_r,
#        txt=paths.rseqc.bamgc_txt
#    benchmark:
#        'benchmark/{sample}_bam_gc.tab'
#    log:
#        'log/{sample}_bam_gc.log'
#    conda:
#        SOURCEDIR+"/../envs/rseqc.yaml"
#    params:
#        sample='{sample}'
#    priority: 1
#    threads: 1
#    shell:
#      '''
#        echo "read_GC.py -i {input.bam} -o rseqc/{params.sample}" | tee {log}
#        read_GC.py -i {input.bam} -o rseqc/{params.sample} 2>> {log}
#
#        ## R script to get txt output info
#        echo "out=as.vector(summary(gc));dta = data.frame('{params.sample}',out[1],out[2],out[3],out[4],out[5],out[6]);write.table(dta,file='{output.txt}',sep="\t",row.names=F,col.names=F,quote=F);" >> {output.r}
#        sed -i "s/pdf/png/g" {output.r} 
#        sed -i 's/main=""/main="{params.sample}"/g' {output.r} 
#        Rscript --vanilla --quiet {output.r}
#
#        ## export rule env details
#        conda env export --no-builds > info/rseqc.info
#      '''
#
