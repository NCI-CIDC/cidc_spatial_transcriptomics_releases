log=file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

source(snakemake@params$functions)

library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(cowplot)

## Resolves error related to future.globals.maxSize.
## Refer to https://github.com/satijalab/seurat/issues/1845
options(future.globals.maxSize = 100000 * 1024^2)
print(getOption("future.globals.maxSize")) 

## Load the filtered Seurat object from the RDS
sample <- readRDS(snakemake@input$rds)
Assays(sample)
sample


## Obtain assay name depending on if sample is SGE or HD
spatial.meta.names <- get_spatial_names(sample, snakemake@params$visium_type, snakemake@params$bin)
assay <- spatial.meta.names$assay
print(assay)

## Perform SCTransform if specified; otherwise, perform default normalization method
if (snakemake@params$sctransform == "TRUE") {
    sample <- SCTransform(sample, assay=assay, verbose=TRUE, conserve.memory=TRUE)
} else {
    sample <- NormalizeData(sample, verbose=TRUE)
    ## Perform imputation with RunALRA() if specified in the config.yaml
    ## This will need further testing to ensure it is being used correclty
    if (snakemake@params$imputation == "TRUE") {
        sample <- RunALRA(sample)
    }
    sample <- FindVariableFeatures(sample, selection.method="vst",
                                   nfeatures=snakemake@params$nfeatures,
                                   verbose=TRUE)
    ## Normalize the data by centering and scaling gene expression values
    sample <- ScaleData(sample)
}

Assays(sample)
sample

## Save normalized data
saveRDS(sample, file=snakemake@output$rds)
