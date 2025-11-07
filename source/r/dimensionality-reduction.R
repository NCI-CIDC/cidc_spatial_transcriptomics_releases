log=file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(Seurat)
library(ggplot2)

## Load the normalized Seurat object from the RDS
sample <- readRDS(snakemake@input$rds)
Assays(sample)
sample

## Run a PCA dimensionality reduction
sample <- RunPCA(sample)

## Scatter plot of the PCA reduction, showing the distribution of samples in the PCA space
plot1 <- DimPlot(sample, reduction="pca", group.by="orig.ident") + ggtitle(NULL)
##  Visualization of the loadings of the first two principal components, indicating the contribution of each feature to these components
plot2 <- VizDimLoadings(sample, dims=1:2)
## Heatmap of the top principal components, displaying the expression of the top principal components across a subset of cells
plot3 <- DimHeatmap(sample, dims=1:9, cells=500, balanced=TRUE, fast=FALSE)
## Generate an elbow plot
plot4 <- ElbowPlot(sample)

## Save plots
ggsave(snakemake@output$pca, plot=plot1)
ggsave(snakemake@output$vizdim, plot=plot2, width=10)
ggsave(snakemake@output$heatmap, plot=plot3, width=13, height=10)
ggsave(snakemake@output$elbow, plot=plot4, bg="white")

## Perform UMAP dimensional reduction on the object using the first 20 principal components
sample <- RunUMAP(sample, dims=1:20)

Assays(sample)
sample

## Save data
saveRDS(sample, file=snakemake@output$rds)
