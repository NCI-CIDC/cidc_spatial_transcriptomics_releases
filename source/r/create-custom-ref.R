log=file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(Seurat)
options(future.globals.maxSize=100000*1024^2)
print(getOption("future.globals.maxSize"))
options(Seurat.object.assay.calcn=TRUE)

create_custom_ref <- function(counts_file_csv, metadata_file_csv) {
    counts <- as.data.frame(data.table::fread(counts_file_csv, sep=",", showProgress=T))
    rownames(counts) <- counts$V1
    counts$V1 <- NULL
    meta.data <- read.table(metadata_file_csv, sep=",")
   
    ## Parse out the annotation levels/categories
    lvls <- colnames(meta.data)
    write(lvls, file=snakemake@output$txt, sep=",")

    ## Generate the reference Seurat object from the counts and metadata csv files
    seurat.obj <- CreateSeuratObject(counts=counts, meta.data=meta.data)
    seurat.obj <- NormalizeData(object=seurat.obj, normalization.method="LogNormalize", scale.factor=10000, margin=1)
    gc()
    seurat.obj <- ScaleData(seurat.obj, features=rownames(seurat.obj), assay="RNA")
    gc()
    seurat.obj <- FindVariableFeatures(object=seurat.obj)
    gc()
    seurat.obj <- RunPCA(seurat.obj, verbose=TRUE)
    gc()
    seurat.obj <- SCTransform(seurat.obj)
    gc()
    seurat.obj <- RunPCA(seurat.obj, verbose=TRUE)
    seurat.obj <- RunUMAP(seurat.obj, dims=1:10, verbose=TRUE)
    gc()

    print(head(seurat.obj@meta.data, 10))
    saveRDS(seurat.obj, file=snakemake@output$rds)
}

create_custom_ref(counts_file_csv=snakemake@input$counts, metadata_file_csv=snakemake@input$metadata)
