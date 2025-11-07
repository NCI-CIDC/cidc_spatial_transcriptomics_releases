log=file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(Seurat)

options(Seurat.object.assay.calcn=TRUE)

## If workflow is set to "10x" then read in the 10X directory provided in the config.yaml.
## Otherwise, the Space Ranger output will be used to generate the Seurat object.

## TODO: Handle case where they provide the directory from Space Ranger
if (snakemake@params$workflow == "10x") {
    sample <- Read10X(snakemake@input$tenx_dir)
} else {
    if (snakemake@params$visium_type == "sge") {
        sample <- Load10X_Spatial(
                                  data.dir=snakemake@params$outs,
                                  filename="filtered_feature_bc_matrix.h5",
                                 )
        cat("Visium Spatial Gene Expression sample detected.\n")
    } else {
        sample <- Load10X_Spatial(
				  data.dir=snakemake@params$outs,
				  filename="filtered_feature_bc_matrix.h5",
				  bin.size=snakemake@params$bin_size
                                 )
        cat("Visium HD Spatial Gene Expression sample detected.\n")
    }
}

sample

## Read in the metadata
metadata.df <- read.csv(snakemake@input$csv, header=TRUE, stringsAsFactors=FALSE)
print(metadata.df)

sample$orig.ident <- snakemake@params$sample

if (nrow(metadata.df) == 0) {
    cat("No additional metadata was identified.\n")
} else {
    ## Select the metadata row for the current sample
    metadata.sample <- subset(metadata.df, samid==snakemake@params$sample)
    metadata.sample <- merge(sample@meta.data, metadata.sample, by.x="orig.ident", by.y="samid", all.x=TRUE)
    ## Add the metadata columns to the Seurat object
    sample <- AddMetaData(sample, metadata=metadata.sample)
    cat("Additional metadata was added to the Seurat object.\n")
}

head(sample@meta.data, 10)
saveRDS(sample, file=snakemake@output$rds)
