log=file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(Seurat)
library(ggplot2)
library(cowplot)

## Load the normalized Seurat object from the RDS
sample <- readRDS(snakemake@input$rds)
Assays(sample)
sample

## Print cell cycle marker genes for S phase and G2/M phase
print(cc.genes.updated.2019$s.genes)
print(cc.genes.updated.2019$g2m.genes)

rownames(sample)
## Perform cell cycle scoring using predefined S phase and G2/M phase gene sets
## Some samples crash at this step, so it is placed in a tryCatch statement
## Refer to https://github.com/satijalab/seurat/issues/1181#issue-415260625
tryCatch({
    sample <- CellCycleScoring(sample, s.features=cc.genes.updated.2019$s.genes,
			       g2m.features=cc.genes.updated.2019$g2m.genes,
			       set.ident=TRUE)
}, error = function(e) {
    message("Error: ", e$message)
    message("An error occurred during cell cycle scoring. Therefore, cell cycle scoring was NOT performed.")
})

head(sample@meta.data)
Assays(sample)
sample

## Save data
saveRDS(sample, file=snakemake@output$rds)
