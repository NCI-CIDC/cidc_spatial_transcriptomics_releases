log=file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

source(snakemake@params$functions)

library(Seurat)
library(dplyr)
library(kableExtra)
library(cowplot)
library(ggplot2)
library(patchwork)
library(stringr)

source(snakemake@params$functions)

sample <- readRDS(snakemake@input$rds)
Assays(sample)
sample

## Set metadata for the sample
metadata <- read.csv(snakemake@input$metadata, header=TRUE, stringsAsFactors=FALSE)
print(metadata)
if (nrow(metadata) > 0) {
    metadata <- metadata %>% filter(samid == snakemake@params$sample)
} else { ## Handle cases when no additional metadata was defined
    metadata <- data.frame(samid=snakemake@params$sample)
}
print(metadata)

## Create run parameters table
parameter <- data.frame(
    Parameter = c("Perform SCTransform normalization",
		  "Remove HLA, IGH, IGK, and IGL genes",
		  "Remove mitochondrial genes",
		  "Remove ribosomal genes",
		  "Remove hemoglobin genes",
		  "Filter spots by minimum nCount threshold",
		  "Filter spots by minimum nFeature threshold",
		  "Annotate with Azimuth reference",
		  "Annotate with Tabula Sapiens reference",
		  "Annotate with custom reference"),
    Value = c(snakemake@params$sctransform,
	      snakemake@params$hla_ig,
	      snakemake@params$mt,
	      snakemake@params$rb,
              snakemake@params$hb,
	      snakemake@params$ncount_min,
	      snakemake@params$nfeature_min,
	      snakemake@params$azimuth,
	      snakemake@params$ts,
	      snakemake@params$custom)
)

## Read in the amount of cells for Seurat object pre- and post- filtering.
nspot.pre <- readLines(snakemake@input$qc_nspot_pre)
print(nspot.pre)
nspot.post <- ncol(sample)
print(nspot.post)

## Heatmap plot
## TODO: Might need to make this an option specifically for the heatmap report section in the config
## Refers to line 794
if (snakemake@params$sctransform == "TRUE") {
    markers <- FindAllMarkers(sample, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25, assay="SCT")
} else {
    markers <- FindAllMarkers(sample, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25, assay="RNA")
}

#print(markers)
top5 <- markers %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)
feat.hm <- top5$gene
#print(top5)
#print(feat.hm)

if (snakemake@params$sctransform == "TRUE") {
    sample <- ScaleData(sample, verbose=FALSE, assay="SCT", features=feat.hm)
    heatmap.plot <- DoHeatmap(sample, features=feat.hm, assay="SCT", slot="scale.data", size=3)
} else {
    sample <- ScaleData(sample, verbose=FALSE, assay="RNA", features=feat.hm)
    heatmap.plot <- DoHeatmap(sample, features=feat.hm, assay="RNA", slot="scale.data", size=3)
}

ggsave(snakemake@output$heatmap, plot=heatmap.plot)

## Azimuth cell type annotation plots
if (snakemake@params$azimuth == "TRUE") {
    ## Sepcify the reference argument used for RunAzimuth()
    azimuth.ref <- snakemake@params$azimuth_ref
    azimuth.lvl <- snakemake@params$azimuth_lvl
    ## Read in levels previously identified in the Azimuth module
    azimuth.lvls <- readLines(snakemake@input$azimuth_txt, n=1)
    azimuth.lvls <- unlist(strsplit(azimuth.lvls, ",\\s*"))
    azimuth.a <- create_dimplot(sample, azimuth.lvls[1])
    azimuth.aa <- create_spatial_dimplot(sample, azimuth.lvls[1], snakemake@params$visium_type)

    ggsave(snakemake@output$azimuth_a, plot=azimuth.a)
    ## width=5, height=4
    ggsave(snakemake@output$azimuth_aa, plot=azimuth.aa)
    azimuth.a <- snakemake@output$azimuth_a
    azimuth.aa <- snakemake@output$azimuth_aa
    ## If there is more than one level, plot the first 2 levels listed
    if (length(azimuth.lvls) > 1) {
        azimuth.b <- create_dimplot(sample, azimuth.lvls[2])
        azimuth.bb <- create_spatial_dimplot(sample, azimuth.lvls[2], snakemake@params$visium_type)
        ggsave(snakemake@params$azimuth_b, plot=azimuth.b)
	ggsave(snakemake@params$azimuth_bb, plot=azimuth.bb)
	azimuth.b <- snakemake@params$azimuth_b
	azimuth.bb <- snakemake@params$azimuth_bb
    } else {
        azimuth.b <- NULL
        azimuth.bb <- NULL
    }
    azimuth.message <- paste0(
        "Cell type annotation was performed using Azimuth with annotation ",
	"level(s): ", azimuth.lvl, ". The reference was set to ",
	azimuth.ref, ". Only the first two annotation levels were used for ",
	"plotting if more than two levels were provided.")
} else {
    azimuth.a <- NULL
    azimuth.aa <- NULL
    azimuth.b <- NULL
    azimuth.bb <- NULL
    azimuth.message <- paste0(
        "Cell type annotation was NOT performed using Azimuth, as the run ",
	"parameters were configured to exclude this step.")
}

## Tabula Sapiens cell type annotation plots
if (snakemake@params$ts == "TRUE") {
    ## Read in the reference type identified in the Tabula Sapiens module
    ts.txt <- readLines(snakemake@input$ts_txt, n=2)
    ts.ref <- ts.txt[-1]
    ts.ref <- gsub("_", " ", ts.ref)
    ts.lvl <- snakemake@params$ts_lvl
    ## Read in categories/levels identified in the Tabula Sapiens module
    ts.lvls <- ts.txt[1]
    ts.lvls <- unlist(strsplit(ts.lvls, ",\\s*"))
    ts.a <- create_dimplot(sample,ts.lvls[1])
    ts.aa <- create_spatial_dimplot(sample, ts.lvls[1], snakemake@params$visium_type)
    ggsave(snakemake@output$ts_a, plot=ts.a)
    ggsave(snakemake@output$ts_aa, plot=ts.aa)
    ts.a <- snakemake@output$ts_a
    ts.aa <- snakemake@output$ts_aa
    ## If there is more than one annotation category/level, plot the first two
    print(length(ts.lvls))
    if (length(ts.lvls) > 1) {
        ts.b <- create_dimplot(sample, ts.lvls[2])
        ts.bb <- create_spatial_dimplot(sample, ts.lvls[2], snakemake@params$visium_type)
        ggsave(snakemake@params$ts_b, plot=ts.b)
	ggsave(snakemake@params$ts_bb, plot=ts.bb)
	ts.b <- snakemake@params$ts_b
	ts.bb <- snakemake@params$ts_bb
    } else {
        ts.b <- NULL
        ts.bb <- NULL
    }
    ts.message <- paste0(
        "Cell type annotation was performed using the Tabula Sapiens ",
	 ts.ref, " reference with annotation categories (or category): ",
	 ts.lvl, ". Only the first two annotation categories were used ",
	 "for plotting if more than two catgories were provided.")
} else {
    ts.a <- NULL
    ts.aa <- NULL
    ts.b <- NULL
    ts.bb <- NULL
    ts.message <- paste0(
        "Cell type annotation was NOT performed using a Tabula Sapiens ",
        "reference, as the run parameters were configured to exclude this step.")
}

## Custom cell type annotation plots
if (snakemake@params$custom == "TRUE") {
    ## Read in the reference type identified in the custom module
    custom.lvls <- readLines(snakemake@input$custom_txt, n=1)
    custom.lvl <- gsub("\\bcustom_", "", custom.lvls)
    custom.lvls <- unlist(strsplit(custom.lvls, ",\\s*"))
    ## Remove "custom_" from the beginning of the annotation category
    custom.a <- create_dimplot(sample, custom.lvls[1])
    custom.aa <- create_spatial_dimplot(sample, custom.lvls[1], snakemake@params$visium_type)
    ggsave(snakemake@output$custom_a, plot=custom.a)
    ggsave(snakemake@output$custom_aa, plot=custom.aa)
    custom.a <- snakemake@output$custom_a
    custom.aa <- snakemake@output$custom_aa
    if (length(custom.lvls) > 1) {
        custom.b <- create_dimplot(sample, custom.lvls[2])
        custom.bb <- create_spatial_dimplot(sample, custom.lvls[2], snakemake@params$visium_type)
	ggsave(snakemake@params$custom_b, plot=custom.b)
	ggsave(snakemake@params$custom_bb, plot=custom.bb)
	custom.b <- snakemake@params$custom_b
	custom.bb <- snakemake@params$custom_bb
    } else {
        custom.b <- NULL
        custom.bb <- NULL
    }
    custom.message <- paste0(
        "Cell type annotation was performed using a custom reference ",
	"with annotation categories (or category): ", custom.lvl, ". ",
	"Only the first two annotation categories were used ",
	"for plotting if more than two catgories were provided.")
} else {
    custom.a <- NULL
    custom.aa <- NULL
    custom.b <- NULL
    custom.bb <- NULL
    custom.message <- paste0(
        "Cell type annotation was NOT performed using a custom reference, ",
	"as the run parameters were configured to exclude this step.")
}

## Spatially variable features
svg <- read.csv(snakemake@input$svg_csv, header=FALSE, stringsAsFactors=FALSE, na.strings="N/A")
print(svg)

if (length(svg) == 1 && is.na(svg[1])) {
    svg.message <- paste0("There were no spatially variable features detected in ",
                          "the tissue section.")
    svg.plot <- NULL
} else {
    svg.message <- paste0("This visualization displays the top spatially variable ",
			  "features detected within the tissue section, ",
			  "highlighting their spatial distribution and ",
			  "expression patterns across the sample.")
    svg.plot <- snakemake@params$svg_plot
}

## Have to set the working directory to the "report" directory or else
## knitr::spin will use the current working directory (PREDIR) and place
## the final report HTML there
setwd(snakemake@params$dir)

## Generate report
test.knit <- file.copy(snakemake@params$template, snakemake@output$r)
knitr::spin(snakemake@output$r, knit=TRUE)
