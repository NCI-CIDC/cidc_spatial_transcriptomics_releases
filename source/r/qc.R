log=file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

source(snakemake@params$functions)

library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)

sample <- readRDS(snakemake@input$rds)
Assays(sample)
sample

head(sample@meta.data, 10)

## Add QC stats for mitochondrial, ribosomal, and hemoglobin gene content
## Genes excluded from the human probe set include ribosomal proteins, so the values will be zero.
## https://www.10xgenomics.com/support/cytassist-spatial-gene-expression/documentation/steps/probe-sets/visium-ffpe-probe-sets-overview
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern="^MT-")
sample[["percent.hb"]] <- PercentageFeatureSet(sample, pattern="^HB[^(P)]")
sample[["percent.rb"]] <- PercentageFeatureSet(sample, pattern="^[Rr][Pp][Ss]|^[Rr][Pp][Ll]")

head(sample@meta.data, 10)

## Obtain the metadata names for nCount and nFeature. The names differ depending on if the sample is SGE or HD
spatial.meta.names <- get_spatial_names(sample, snakemake@params$visium_type, snakemake@params$bin_size)
ncount <- spatial.meta.names$ncount
nfeature <- spatial.meta.names$nfeature

## Print out the entries that have 0 as nCount
ncount.zero <- sample@meta.data[[ncount]] == 0
if (any(ncount.zero)) {
    subset_data <- sample[, ncount.zero]
    print(subset_data@meta.data)
} else {
    cat("There are no entries in the Seurat object where nCount equals 0.\n")
}

## Handle cases when NaN values get generated from PercentageFeatureSet because nCount equals 0
sample@meta.data$percent.mt[which(is.na(sample$percent.mt))] = 0
sample@meta.data$percent.hb[which(is.na(sample$percent.hb))] = 0
sample@meta.data$percent.rb[which(is.na(sample$percent.rb))] = 0

## Generate summary QC plot
create_qc_plot(sample, ncount, nfeature, "percent.mt", "percent.hb", "percent.rb",
	       snakemake@output$plot, "coral2")

## Generate spatial plots for each individual feature
create_spatial_plot(ncount, snakemake@output$ncount, "coral2", snakemake@params$visium_type)
create_spatial_plot(nfeature, snakemake@output$nfeature, "coral2", snakemake@params$visium_type)
create_spatial_plot("percent.mt", snakemake@output$mt, "coral2", snakemake@params$visium_type)
create_spatial_plot("percent.hb", snakemake@output$hb, "coral2", snakemake@params$visium_type)
create_spatial_plot("percent.rb", snakemake@output$rb, "coral2", snakemake@params$visium_type)

## Write out spot count of Seurat object for individual report
nspot <- length(rownames(sample@meta.data))
write.table(nspot, file=snakemake@output$nspot, row.names=FALSE, col.names=FALSE)

## Save Seurat object
saveRDS(sample, file=snakemake@output$rds)
