log=file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(Seurat)

## Identify the total amount of samples to process 
total <- length(snakemake@input) - 1
print(total)

sample.list <- list()

## Load the Seurat objects and put them in a list
for (i in 1:total) {
    sample.name <- paste0("sample.", i)
    assign(sample.name, readRDS(snakemake@input[[i]]))
    sample.list[[i]] <- get(sample.name)
}

print(sample.list)

## Rename certain parts of the Seurat objects, so they are consistent
for (i in seq_along(sample.list)) {
    obj <- sample.list[[i]]
    obj@images$slice1
    assay.name <- names(obj@assays)[grep("Spatial", names(obj@assays))]
    DefaultAssay(obj) <- assay.name
    if (names(obj@assays)[grep("Spatial", names(obj@assays))] != "Spatial") {
        assay.name <- names(obj@assays)[grep("Spatial", names(obj@assays))]
        print(assay.name)
	assay.suffix <- paste0(".", strsplit(assay.name, split = ".", fixed=T)[[1]][2])
	print(assay.suffix)
	obj@assays$Spatial <- obj@assays[[assay.name]]
	obj@assays[assay.name] <- NULL
	obj@images$slice1 <- obj@images[[paste0("slice1", assay.suffix)]]
	DefaultAssay(obj) <- "Spatial"
	obj@images[paste0("slice1", assay.suffix)] <- NULL
	obj@images$slice1@assay <- "Spatial"
    }
    obj@images[[unique(obj$orig.ident)]] <- obj@images$slice1
    obj@images$slice1 <- NULL
    sample.list[[i]] <- obj
}

merged.obj <- Reduce(function(x, y) merge(x, y), sample.list)
merged.obj <- JoinLayers(merged.obj, assay="Spatial")
DefaultAssay(merged.obj) <- "Spatial"

merged.obj

## Normalization
merged.obj <- NormalizeData(merged.obj)
merged.obj <- ScaleData(merged.obj, features=rownames(merged.obj))
merged.obj <- FindVariableFeatures(merged.obj, selection.method="vst")

## Dimensionality Reduction
merged.obj <- RunPCA(merged.obj, features=VariableFeatures(object=merged.obj))
merged.obj <- RunUMAP(merged.obj, dims=1:10)

## Clustering
merged.obj <- FindNeighbors(merged.obj, dims=1:10)
merged.obj <- FindClusters(merged.obj, resolution=0.01, set.seed(length(merged.obj)))
merged.obj <- FindClusters(merged.obj, resolution=0.05, set.seed(length(merged.obj)))
merged.obj <- FindClusters(merged.obj, resolution=0.1, set.seed(length(merged.obj)))
merged.obj <- FindClusters(merged.obj, resolution=0.2, set.seed(length(merged.obj)))
merged.obj <- FindClusters(merged.obj, resolution=0.4, set.seed(length(merged.obj)))
merged.obj <- FindClusters(merged.obj, resolution=0.8, set.seed(length(merged.obj)))
merged.obj <- FindClusters(merged.obj, resolution=1.2, set.seed(length(merged.obj)))
merged.obj <- FindClusters(merged.obj, resolution=1.8, set.seed(length(merged.obj)))
merged.obj <- FindClusters(merged.obj, resolution=2.4, set.seed(length(merged.obj)))
merged.obj <- FindClusters(merged.obj, resolution=3.0, set.seed(length(merged.obj)))

merged.obj <- RunTSNE(merged.obj, assay="Spatial", reduction="pca", dims=1:10)
merged.obj <- JoinLayers(merged.obj, assay="Spatial")

merged.obj

## Save data
saveRDS(merged.obj, file=snakemake@output$rds)
