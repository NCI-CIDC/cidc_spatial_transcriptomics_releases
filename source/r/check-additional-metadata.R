log=file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Read in the metadata from the configuration file sample_metadata.csv
metadata.df <- read.csv(snakemake@input[[1]], header=TRUE, stringsAsFactors=FALSE, comment.char="#")

print(metadata.df)

# Drop columns that are not needed for the Seurat object meta.data
metadata.df <- subset(metadata.df, select = -c(fastqs_dir, cytaimage, loupe_alignment,
					       slide, area, probe_set_version, dapi_index,
					       image, darkimage, colorizedimage))
print(metadata.df)
print(ncol(metadata.df))

if (ncol(metadata.df) > 1) {
    cat("Additional metadata was identified in the sample_metadata.csv.\n")
    write.csv(metadata.df, snakemake@output[[1]], row.names=FALSE)
} else {
    cat("No additional metadata was identified in the sample_metadata.csv.\n")
    writeLines("N/A", snakemake@output[[1]])
}

## Remaining column names are split by data
## TODO: Output those column names in a different document


## TODO: Logic for when there is no split.by data
