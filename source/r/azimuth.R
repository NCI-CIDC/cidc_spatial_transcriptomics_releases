log=file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)
library(ggplot2)

options(future.globals.maxSize = 8000 * 1024^2)

## Load the normalized Seurat object from the RDS
sample <- readRDS(snakemake@input$rds)
Assays(sample)
sample

## Run Azimuth annotation with specified reference
sample <- RunAzimuth(sample, reference=snakemake@params$ref)

head(sample@meta.data, 10)

## Parse out the reference type
type <- sub("ref$", "", snakemake@params$ref)

lvls <- unlist(strsplit(snakemake@params$lvl, ", "))
#lvls <-  unlist(strsplit("celltype.l1", ", "))

## Vector to hold the first two annotation levels listed
updated.cell.names <- c()

## Rename annotation columns and add a score bins column for each annotation level
for (i in seq_along(lvls)) {
    lvl <- lvls[i]
    cell <- paste0("predicted.", lvl)
    score <- paste0("predicted.", lvl, ".score")

    ## Rename the annotation columns to include the reference type and the sample type
    cell.rename <- paste("Azimuth", type, lvl, sep=".")
    score.rename <- paste("Azimuth", type, lvl, "score", sep=".")
    names(sample@meta.data)[names(sample@meta.data) == cell] <- cell.rename
    names(sample@meta.data)[names(sample@meta.data) == score] <- score.rename

    ## Save the updated column name for the first two annotation levels in the list
    if (i <= 2) {
        updated.cell.names <- c(updated.cell.names, cell.rename)
    }

    ## Add score bins as a column
    bin <- paste("Azimuth", type, lvl, "score.bins", sep =".")
    sample@meta.data[[bin]] <- cut(round(sample@meta.data[[score.rename]], 2), breaks=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),
                                        labels=c("0-10","11-20","21-30","31-40","41-50","51-60","61-70","71-80","81-90","91-100"))
}

head(sample@meta.data, 10)

print(length(lvls))

## Save the renamed levels for use in the report
updated.cell.names.list <- paste(updated.cell.names, collapse=", ")
print(updated.cell.names.list)
writeLines(updated.cell.names.list, snakemake@output$txt)

sample

## Save data
saveRDS(sample, file=snakemake@output$rds)
