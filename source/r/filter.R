log=file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

source(snakemake@params$functions)

library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)

## Filter based on min.cells and min.features 
sample <- readRDS(snakemake@input[[1]])
Assays(sample)
sample

all.genes <- rownames(sample)
remove.genes <- character(0)
## Remove HLA, IGK, IGH, and IGL genes from Seurat object
if (snakemake@params$hla_ig == "TRUE") {
    remove.genes <- c(remove.genes, all.genes[grepl("^[Ii][Gg][Hh]", all.genes) &
		      !grepl("^IGHMBP", all.genes) |
	              grepl("^[Ii][Gg][Kk]", all.genes) |
		      grepl("^[Ii][Gg][Ll]", all.genes) |
		      grepl("^[Hh][Ll][Aa]", all.genes)])
}

## Remove mitochondrial genes from Seurat object
if (snakemake@params$mt == "TRUE") {
    remove.genes <- c(remove.genes, all.genes[grepl("^MT-", all.genes)])
}

## Remove ribosomal genes from Seurat oject
if (snakemake@params$rb == "TRUE") {
    remove.genes <-c(remove.genes, all.genes[grepl("^[Rr][Pp][Ss]", all.genes) |
		     grepl("^[Rr][Pp][Ll]", all.genes)])
}

## Remove hemoglobin genes from Seurat object
if (snakemake@params$hb == "TRUE") {
    remove.genes <- c(remove.genes, all.genes[grepl("^HB[^(P)]", all.genes)])
}

print(all.genes)
print(remove.genes)
if (length(remove.genes) > 0) {
    sample <- sample[all.genes[!all.genes %in% remove.genes],]
}

sample

## Specify names of nCount and nFeature depending on if it is an SGE or HD sample
spatial.meta.names <- get_spatial_names(sample, snakemake@params$visium_type, snakemake@params$bin_size)
ncount <- spatial.meta.names$ncount
nfeature <- spatial.meta.names$nfeature

## Obtain QC stats for mitochondrial, ribosomal, and hemoglobin genes after filtering
sample[["percent.mt.filtered"]] <- PercentageFeatureSet(sample, pattern="^MT-")
sample[["percent.rb.filtered"]] <- PercentageFeatureSet(sample, pattern="^[Rr][Pp][Ss]|^[Rr][Pp][Ll]")
sample[["percent.hb.filtered"]] <- PercentageFeatureSet(sample, pattern="^HB[^(P)]")

## Subset Seurat object based on nCount and nFeature minimums
keep.spots <- (sample@meta.data[[ncount]] > snakemake@params$ncount_min) &
	      (sample@meta.data[[nfeature]] > snakemake@params$nfeature_min)
sample <- sample[, keep.spots]

## Generate summary QC plot
create_qc_plot(sample, ncount, nfeature, "percent.mt.filtered", "percent.hb.filtered",
	       "percent.rb.filtered", snakemake@output$plot, "cornflowerblue")

## Generate spatial plots for each individual feature
create_spatial_plot(ncount, snakemake@output$ncount, "cornflowerblue", snakemake@params$visium_type)
create_spatial_plot(nfeature, snakemake@output$nfeature, "cornflowerblue", snakemake@params$visium_type)
create_spatial_plot("percent.mt.filtered", snakemake@output$mt, "cornflowerblue", snakemake@params$visium_type)
create_spatial_plot("percent.hb.filtered", snakemake@output$hb, "cornflowerblue", snakemake@params$visium_type)
create_spatial_plot("percent.rb.filtered", snakemake@output$rb, "cornflowerblue", snakemake@params$visium_type)

## Write out spot count of Seurat object for individual report
nspot <- length(rownames(sample@meta.data))
write.table(nspot, file=snakemake@output$nspot, row.names=FALSE, col.names=FALSE)

head(sample@meta.data, 10)

## Save filtered data
saveRDS(sample, file=snakemake@output$rds)
