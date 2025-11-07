log=file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(remotes)
library(httr)

## Update packages with versions that are not available from conda yet
packages <- c(
	      "dqrng"="0.4.1",
	      "abind"="1.4-8",
              "openssl"="2.3.1",
	      "curl"="6.1.0",
	      "data.table"="1.16.4",
	      "Rnanoflann"="0.0.3",
	      "Rfast"="2.1.5.1",
	      "Rfast2"="0.1.5.2"
	     )

## Loop through each package and version and install to the R environment
for (pkg in names(packages)) {
    version <- packages[pkg]
    cat("Installing", pkg, "version", version, "\n")
    ## Install the package with the specified version
    devtools::install_version(pkg, version=version, repos="http://cran.rstudio.com/", force=TRUE, upgrade="never")
}

## API rate limit for unauthorized requests can be exceeded resulting in error 403.
## This error may occur if the script is executed too many times consecutively.
## Monitor the remaining API requests by checking the "remaining" value.
## The rate limit is consumed with each call to remotes::install_github.
## Refer to https://github.com/r-lib/remotes/issues/210 for example of error.
response <- GET("https://api.github.com/rate_limit")
print(response)

## Install SeuratWrappers (https://github.com/satijalab/seurat-wrappers)
remotes::install_github("satijalab/seurat-wrappers@8d46d6c", force=TRUE, upgrade="never")

response <- GET("https://api.github.com/rate_limit")
print(response)

## Install version 0.2.2.9001 of SeuratData to resolve issue with version 2.1.0
## Refer to https://github.com/satijalab/seurat/issues/6434 and
## https://github.com/satijalab/azimuth/issues/121
packageVersion("SeuratData")
devtools::install_github("satijalab/seurat-data@4dc08e0", force=TRUE, upgrade="never")
packageVersion("SeuratData")

## Output the installed package versions
installed.packages <- as.data.frame(installed.packages()[, c("Package", "Version")])
print(installed.packages)
write.table(installed.packages, file=snakemake@output$csv, row.names=FALSE, quote=FALSE, sep="\t")
