#' ---
#' title: "Spatial Transcriptomics `r snakemake@params$sample` Report"
#' ---
#'
#' ## Sample Metadata
#' The table contains detailed metadata for the sample as defined in the sample_metadata.csv file.
#' ```{r, echo=FALSE}
#' kbl(metadata) %>%
#'   kable_styling() %>%
#'   scroll_box(width = "100%")
#' ```
#'
#' ## Run Parameters
#' ### Run Configuration Overview
#' The following processing steps and options can be defined by configuring the corresponding settings in the config.yaml file.
#'
#' ```{r, echo=FALSE}
#' kbl(parameter)
#' ```
#'
#' ## QC Summary
#' ### Pre-filtering
#' The plotting is performed before any filtering steps are applied to the sample. The sample contained `r nspot.pre` spots, and this number reflects the total spot count prior to any quality control and gene filtering.
#' ![](`r snakemake@input$qc_plot_pre`)
#'
#' The plot visualizes the distribution and spatial organization of gene expression counts (nCount) across the tissue section.
#' ![](`r snakemake@input$ncount_pre`)
#'
#' The plot visualizes the spatial distribution of mitochondrial gene expression across the tissue section.
#' ![](`r snakemake@input$mt_pre`)
#'
#' ### Post-filtering
#' The plotting is performed after filtering steps are applied to the sample. The sample contains `r nspot.post` cells, and this number reflects the total spot count after quality control and gene filtering.
#' ![](`r snakemake@input$qc_plot_post`)
#'
#' The plot visualizes the distribution and spatial organization of gene expression counts (nCount) across the tissue section.
#' ![](`r snakemake@input$ncount_post`)
#'
#' The plot visualizes the spatial distribution of mitochondrial gene expression across the tissue section.
#' ![](`r snakemake@input$mt_post`)
#'
#' ## Cell Type Annotation
#' ### Azimuth
#' Azimuth is an automated tool for cell type annotation that utilizes a pre-annotated reference single-cell dataset and a feature-barcode matrix from the query dataset as inputs. It assigns cell types by comparing the gene expression profiles of individual cells in the query dataset to those in the reference dataset.
#' Visit the Azimuth reference page for a list of references available: https://azimuth.hubmapconsortium.org/references/
#'
#' `r azimuth.message`
#' `r if (!is.null(azimuth.a)) paste0("![](", azimuth.a, ")")`
#' `r if (!is.null(azimuth.a)) paste0("![](", azimuth.aa, ")")`
#' `r if (!is.null(azimuth.b)) paste0("![](", azimuth.b, ")")`
#' `r if (!is.null(azimuth.b)) paste0("![](", azimuth.bb, ")")`
#'
#' ### Tabula Sapiens
#' Tabula Sapiens is a comprehensive single-cell reference atlas that provides detailed annotations of cell types across multiple human tissues. A reference is generated from the Tabula Sapiens atlas and used for cell type annotation. For more information, visit the Tabula Sapiens resource page: https://tabula-sapiens.sf.czbiohub.org/
#'
#' `r ts.message`
#' `r if (!is.null(ts.a)) paste0("![](", ts.a, ")")`
#' `r if (!is.null(ts.aa)) paste0("![](", ts.aa, ")")`
#' `r if (!is.null(ts.b)) paste0("![](", ts.b, ")")`
#' `r if (!is.null(ts.bb)) paste0("![](", ts.bb, ")")`
#'
#' ### Custom
#' A custom reference is generated from provided counts.csv and metadata.csv files. This reference is then used to perform cell type annotation.
#'
#' `r custom.message`
#' `r if (!is.null(custom.a)) paste0("![](", custom.a, ")")`
#' `r if (!is.null(custom.aa)) paste0("![](", custom.aa, ")")`
#' `r if (!is.null(custom.b)) paste0("![](", custom.b, ")")`
#' `r if (!is.null(custom.bb)) paste0("![](", custom.bb, ")")`
#'
#' ## Clustering
#' The plot is a UMAP visualization of the data where each point represents a cell, and cells are grouped into clusters based on similar gene expression profiles. The clusters, labeled numerically, highlight distinct cell types or states.
#' ![](`r snakemake@input$cluster_umap`)
#' The plotting visualizes the spatial location of different clusters individually with a separate plot for each cluster. This allows for a clear examination of the spatial location and organization of each cluster within the tissue.
#' ![](`r snakemake@input$cluster_individual`)
#'
#' ## Heatmap Cluster Genes
#' The plot is a heatmap visualization of the top marker genes identified for each cluster in the data.
#' ![](`r snakemake@output$heatmap`)
#'
#' The table presents the top 5 genes identified for each cluster along with their associated statistical and expression metrics.
#' ```{r, echo=FALSE}
#' kbl(cbind(top5)) %>%
#'     kable_styling() %>%
#'     scroll_box(height = "250px")
#' ```
#'
#' ## Spatially Variable Features
#' `r svg.message`
#' `r if (!is.null(svg.plot)) paste0("![](", svg.plot, ")")`
