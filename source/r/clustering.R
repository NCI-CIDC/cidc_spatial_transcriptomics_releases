log=file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(Seurat)
library(ggplot2)
library(patchwork)

## Load the Seurat object from the RDS
sample <- readRDS(snakemake@input$rds)
sample

## Set the seeed as the number of samples (SAMID) to process for the batch
seed <- length(snakemake@params$samid)
print(snakemake@params$samid)
print(seed)

## Compute the k.param nearest neighbors
sample <- FindNeighbors(sample, dims=1:10)

## Identify clusters of cells by a shared nearest neighbor (SNN) modularity
## optimization based clustering algorithm
sample <- FindClusters(sample, set.seed(seed), resolution=0.25)

## Generate clustering plots for visualization
clusters <- levels(sample$seurat_clusters)
print(clusters)
plot1 <- DimPlot(sample, group.by="seurat_clusters", label=TRUE,
		 raster=FALSE, repel=TRUE) + labs(title=NULL) +
                 theme(aspect.ratio=1)
if (snakemake@params$visium_type == "hd") {
    point.size.factor <- 10
} else {
    point.size.factor <- 1.6
}
plot2 <- SpatialDimPlot(sample, label=FALSE,
			pt.size.factor=point.size.factor) +
                        theme(legend.title=element_blank(),
			legend.key.spacing.y=unit(0, "cm")) +
	                theme(aspect.ratio=1)
plot1 <- wrap_plots(plot1, plot2)

cells <- CellsByIdentities(sample, idents=clusters)
plot3 <- SpatialDimPlot(sample,
		        cells.highlight=cells[setdiff(names(cells), "NA")],
		        cols.highlight=c("#FFFF00", "grey50"),
		        facet.highlight=T, combine=T) + NoLegend()
ggsave(snakemake@output$umap, plot=plot1, dpi=300, width=12, height=6)
ggsave(snakemake@output$individual, plot=plot3)

head(sample@meta.data, 10)
sample

## Save data
saveRDS(sample, file=snakemake@output$rds)
