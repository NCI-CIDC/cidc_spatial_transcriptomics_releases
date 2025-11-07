log=file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(Seurat)
library(ggplot2)
library(patchwork)

options(width = 200)

## Load the Seurat object from the RDS
sample <- readRDS(snakemake@input$rds)
sample

if (snakemake@params$sctransform == "TRUE") {
    assay.type <- "SCT"
} else {
    assay.type <- "RNA"
}

## Identify features whose variability in expression can be explained
## to some degree by spatial location
sample <- FindSpatiallyVariableFeatures(sample, assay=assay.type,
					features=VariableFeatures(sample),
					selection.method="moransi", verbose=TRUE,
					x.cuts=50, y.cuts=50)

## TODO: Code that handles when spatially variable features are not identified needs to be tested
## and when the amount of variable features is less than 10

## Filter Moran's I feature table
meta <- sample@assays[[assay.type]]@meta.features
rank <- meta[!is.na(meta$moransi.spatially.variable.rank), ] ## Remove NA rows
rank <- rank[order(rank$moransi.spatially.variable.rank), ] ## Sort descending
rank ## Print out the full table

## Subset for the features that were identified as spatially variable.
## 2000 is the default for the number of spatially variable features that
## should be identified.
true.features <- rownames(rank)[rank$moransi.spatially.variable == TRUE]
length(true.features)

if (length(true.features) >= 10) {
    nfeatures <- 10
} else {
    nfeatures <- length(true.features)
}

print(nfeatures)

if (nfeatures > 0) {
    top.features <- true.features[1:nfeatures]
    print(top.features)

    ## Generate an all-encompassing PNG that includes each feature plot
    all.feature.plot <- SpatialFeaturePlot(sample,features=top.features,
					   pt.size.factor=1.6,
				           ncol=3, alpha=c(0.9, 1), image.alpha=0.5)
    ## TODO: Code might need to be dynamic for the plotting width and height depending
    ## on how many features end up getting plotted
    ggsave(snakemake@params$png, plot=all.feature.plot, dpi=300, width=10, height=15)

    ## Generate an individual SpatialFeaturePlot for each feature
    for (feature in top.features) {
        feature.plot <- SpatialFeaturePlot(sample, features=feature,
					   pt.size.factor=1.6,
					   ncol=1, alpha=c(0.9, 1), image.alpha=0.5)
	plot.path <- paste0(snakemake@params$path, snakemake@params$sample,
			    "_spatial_variable_feature_", feature, ".png")
	print(plot.path)
        ggsave(plot.path, dpi=300, width=6, height=6)
    }
    
    ## Write top features to a CSV file that will be used for the sample report
    top.features <- paste(top.features, collapse=",")
    print(top.features)
    write.table(top.features, snakemake@output$csv, row.names=FALSE, col.names=FALSE, quote=FALSE)

} else { ## Handle cases when no spatially variable features are identified
    cat("No spatially variable features were identified.\n")
    writeLines("N/A", snakemake@output$csv)
}

head(sample@meta.data, 10)
sample

## Save data
saveRDS(sample, file=snakemake@output$rds)
