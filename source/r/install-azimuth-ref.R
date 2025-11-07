log=file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(SeuratData)

## Print the available references
available_data <- AvailableData()
print(available_data[grep("Azimuth", available_data[, 3]), 1:3])

## Install reference
InstallData(snakemake@params[[1]], force.reinstall=TRUE)

file.create(snakemake@output[[1]])
