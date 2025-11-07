log=file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(Seurat)
library(patchwork)
library(ggplot2)

options(future.globals.maxSize = 8000 * 1024^2)

## Load the normalized Seurat object from the RDS
ref <- readRDS(snakemake@input$ref)
sample <- readRDS(snakemake@input$rds)

## Set the normalization method used
if (snakemake@params$sctransform == 'TRUE') {
    norm.method <- "SCT"
} else {
    DefaultAssay(ref) <- "RNA"
    norm.method <- "LogNormalize"
}

Assays(sample)
sample
Assays(ref)
ref

head(sample@meta.data, 10)
head(ref@meta.data, 10)

## Find a set of anchors between a reference and query object
anchors <- FindTransferAnchors(reference=ref, query=sample, normalization.method=norm.method,
			       reference.reduction="pca", dims=1:50, k.anchor=snakemake@params$k_anchor)

## Annotation categories/levels for Tabula Sapiens reference
lvls <- unlist(strsplit(snakemake@params$lvl, ", "))

refdata <- setNames(as.list(lvls), paste0("TS_", lvls))
print(refdata)

sample <-MapQuery(anchorset=anchors, query=sample, reference=ref,
		  refdata=refdata, reference.reduction="pca")

## Parse out the reference type from the path 
## The reference file name must be in the format: reference_SCT_TS_{type}.rds
type <- sub(".*TS_(.*)\\.rds$", "\\1", snakemake@input$ref)

## Vector to hold the first two annotation levels listed
updated.cell.names <- c()

head(sample@meta.data, 10)

## Rename annotation columns and add score bins column
for (i in seq_along(lvls)) {
    lvl <- lvls[i]
    cell <- paste0("predicted.TS_", lvl)
    score <- paste0("predicted.TS_", lvl, ".score")

    ## Rename the annotation columns to include the sample type
    cell.rename <- paste("TS", lvl, type, sep="_")
    score.rename <- paste0("TS_", lvl, "_", type, ".score")
    names(sample@meta.data)[names(sample@meta.data) == cell] <- cell.rename
    names(sample@meta.data)[names(sample@meta.data) == score] <- score.rename

    ## Save the updated column name for the first two annotation levels in the list
    if (i <= 2) {
        updated.cell.names <- c(updated.cell.names, cell.rename)
    }

    ## Add score bins as a column
    bin <- paste0(score.rename, ".bins")
    sample@meta.data[[bin]] <- cut(round(sample@meta.data[[score.rename]], 2), breaks=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),
				   labels=c("0-10","11-20","21-30","31-40","41-50","51-60","61-70","71-80","81-90","91-100"))
}

## Save the annotation categories/levels and the reference type
updated.cell.names.list <- paste(updated.cell.names, collapse=", ")
print(updated.cell.names.list)
print(type)
writeLines(c(updated.cell.names.list, type), snakemake@output$txt)

head(sample@meta.data, 10)
sample

## Save data
saveRDS(sample, file=snakemake@output$rds)
