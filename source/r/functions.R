## Specify names of nCount, nFeature, and assay depending on if it is an SGE or HD sample
get_spatial_names <- function(sample, visium.type, bin.size) {
    if (visium.type == "sge") {
        ncount <- "nCount_Spatial"
        nfeature <- "nFeature_Spatial"
	assay <- "Spatial"
    } else {
	suffix <- sprintf("%03dum", bin.size)
        ncount <- paste0("nCount_Spatial.", suffix)
	nfeature <- paste0("nFeature_Spatial.", suffix)
	assay <- paste0("Spatial.", suffix)
    }
    return(list(ncount=ncount, nfeature=nfeature, assay=assay))
}

## Generate summary QC plot
create_qc_plot <- function(sample, nfeature, ncount, mt, hb, rb, plot.path, color) {
    plot1 <- VlnPlot(sample, features=c(nfeature, ncount, mt, hb, rb),
	          ncol=5, group.by="orig.ident", pt.size=0, cols=color) &
                  theme(plot.title=element_text(size=8), axis.text.x=element_text(size=8),
		  axis.title.x = element_text(size=9))
    plot2 <- FeatureScatter(sample, feature1=ncount, feature2=nfeature,
			    group.by="orig.ident", pt.size=0.5, cols=color) +
	                 theme(legend.title=element_text(size=9), legend.text=element_text(size=8))
    plot3 <- ggdraw() + draw_plot(plot1, x=0, y=1/2, width=1, height=1/2) +
	  draw_plot(plot2, x=0, y=0, width=1, height=1/2)
    ggsave(plot.path, plot=plot3, bg="white")
}

## Generate a violin plot and a feature plot from spatial data
create_spatial_plot <- function(feature, plot.path, color, visium.type) {
    if (visium.type == "hd") {
        point.size.factor <- 1.6
    } else {
        point.size.factor <- 1.6
    }
    plot1 <- VlnPlot(sample, features=feature, group.by="orig.ident", pt.size=0, cols=color) +
                     NoLegend() + theme(aspect.ratio=1)
    plot2 <- SpatialFeaturePlot(sample, features=feature, pt.size.factor=point.size.factor) +
	                        theme(legend.position="right") +
                                theme(aspect.ratio=1)
    plot3 <- wrap_plots(plot1, plot2)
    ggsave(plot.path, plot=plot3, dpi=300, width=12, height=6)
}

## Construct the specified column name to be used
column_construct <- function(ref.type, azimuth, ts, lvl) {
    if (ref.type == "azimuth") {
        ref <- sub("ref$", "", azimuth)
        column <- paste("Azimuth", ref, lvl, sep=".")
    } else if (ref.type == "ts") {
        ref <- sub(".*TS_(.*)\\.rds$", "\\1", ts)
	column <- paste("TS", lvl, ref, sep="_")
    } else if (ref.type == "custom") {
        column <- paste0("custom_", lvl)
    }
    return (column)
}

## Check to see if column exists in Seurat object
column_exists <- function(obj, col, file=NULL) {
    col.check <- col %in% colnames(obj@meta.data)
    if (is.null(file)) {
        print(sprintf("The column '%s' was found in the Seurat object: %s", col, col.check))
    } else {
        cat(sprintf("\nThe column \"%s\" was found in the Seurat object: %s\n", col, col.check), file=file)
    }
    return (col.check)
}

## Return proper name of reference
ref_return <- function(ref.type) {
    if (ref.type == "azimuth") {
        ref <- "Azimuth"
    } else if (ref.type == "ts") {
        ref <- "Tabula Sapiens"
    } else if (ref.type == "custom") {
        ref <- "Custom"
    }
    return (ref)
}

## Select metadata names/columns in the Seurat object
select_metadata_names <- function(obj) {
    print(colnames(obj@meta.data))
    selected.metadata <- colnames(obj@meta.data)
    count.unique <- rapply(obj@meta.data, function(x) length(unique(x))) < 31
    selected.metadata <- selected.metadata[count.unique]
    selected.metadata <- selected.metadata[!grepl("score\\.bins", selected.metadata)]
    selected.metadata <- selected.metadata[!grepl("pANN_", selected.metadata)]
    selected.metadata <- selected.metadata[!grepl("\\.Score", selected.metadata)]
    selected.metadata <- selected.metadata[!grepl("\\.score", selected.metadata)]
    selected.metadata <- selected.metadata[!grepl("percent.", selected.metadata)]
    selected.metadata <- selected.metadata[!grepl("nFeature_", selected.metadata)]
    selected.metadata <- selected.metadata[!grepl("nCount_", selected.metadata)]
    print(selected.metadata)
    return(selected.metadata)
}

## Creates a density plot of annotation scores with a vertical red line indicating the filtering cutoff,
## allowing for a visual comparison of the score distribution against the threshold
density_plot <- function(obj, col.score, cutoff) {
    plot <- tibble(value = obj@meta.data[,col.score]) %>%
    ggplot(aes(value)) +
    geom_density(fill = "grey", alpha = 0.8) +
    geom_vline(xintercept = cutoff, color = "red") +
    scale_x_continuous(name = col.score) +
    scale_y_continuous(name = "Density") +
    theme_bw()
    return (plot)
}

## Creates percentage plot
percentage_plot <- function(seurat.obj, group.by, split.by=NA,
			   text.angle=45, text.size=5, file.name="plot.png",
			   width=16, height=10, title.plot="",
                           order.by=NA, plot.file="png") {
    ## Generating custom color palettes for the plot
    custom_colors <- list()
    colors_dutch <- c(
        "#FFC312","#C4E538","#12CBC4","#FDA7DF","#ED4C67",
        "#F79F1F","#A3CB38","#1289A7","#D980FA","#B53471",
        "#EE5A24","#009432","#0652DD","#9980FA","#833471",
        "#EA2027","#006266","#1B1464","#5758BB","#6F1E51"
    )
    colors_spanish <- c(
        "#40407a","#706fd3","#f7f1e3","#34ace0","#33d9b2",
	"#2c2c54","#474787","#aaa69d","#227093","#218c74",
	"#ff5252","#ff793f","#d1ccc0","#ffb142","#ffda79",
	"#b33939","#cd6133","#84817a","#cc8e35","#ccae62"
    )
    custom_colors$discrete <- c(colors_dutch, colors_spanish)

    ## TODO: Check if custom_colors$cell_cyle  actually get used. Can potentially be removed.
    custom_colors$cell_cycle <- setNames(
        c("#45aaf2", "#f1c40f", "#e74c3c", "#7f8c8d"),
	c("G1", "S", "G2M", "-")
    )

    ## Extract metadata from the Seurat object
    meta.data <- seurat.obj@meta.data
    ## If no specific split-by group is provided, create a single sample group (dataset)
    if(is.na(split.by[1])) {
        meta.data$sample <- "dataset"
        meta.data$sample <- factor(meta.data$sample, levels="dataset")
        meta.data$data.plot <- seurat.obj@meta.data[, group.by]
	#print(meta.data$data.plot)
    } else {
	## If split.by is specified, set sample groups accordingly
        if(is.na(order.by[1])) {
            meta.data$sample <- factor(seurat.obj@meta.data[, split.by], levels=unique(seurat.obj@meta.data[, split.by]))
        } else { ## TODO: Check if this code is actually implemented and needed
            meta.data$sample <- factor(seurat.obj@meta.data[, split.by], levels=order.by)
        }
        ## Extract data for plotting based on the group.by variable
        meta.data$data.plot <- seurat.obj@meta.data[, group.by]
    }

    ## Group cells by sample and cell type then calculate cell counts for each
    table.samples.by.cell.type <- meta.data %>%
        dplyr::group_by(sample, data.plot) %>%
        dplyr::summarize(count=n()) %>%
        tidyr::spread(data.plot, count, fill = 0) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(total_cell_count=rowSums(.[c(2:ncol(.))])) %>%
        dplyr::select(c(sample, "total_cell_count", dplyr::everything()))

    ## Count the number of cells in each sample group for labeling in the plot
    temp.labels <- meta.data %>% group_by(sample) %>% tally()

    ## Prepare the data for plotting as a stacked bar chart
    p1_l1 <- table.samples.by.cell.type %>%
        dplyr::select(-c("total_cell_count")) %>%
	reshape2::melt(id.vars="sample") %>%
	mutate(sample=factor(sample, levels=levels(meta.data$sample))) %>%
        ggplot(aes(sample, value)) +
	geom_bar(aes(fill=variable), position="fill", stat="identity") +
	geom_text(data=temp.labels, aes(x=sample, y=Inf,
            label=paste0("n = ", format(n, big.mark = ",", trim=TRUE)), vjust=-1),
            color="black", size=4) +
        scale_fill_manual(name="Cell type", values=custom_colors$discrete) +
	scale_y_continuous(name="Percentage [%]", labels = scales::percent_format(scale=1),
	    expand = c(0.01, 0)) +
        coord_cartesian(clip="off") +
	theme_bw() +
	theme(legend.position="right",
	    plot.title=element_text(hjust=0.5),
	    text=element_text(size=16),
	    panel.grid.major=element_blank(),
	    panel.grid.minor=element_blank(),
	    axis.title.x=element_blank(),
	    axis.text.x=element_text(angle=text.angle, hjust=1, vjust=1, size=text.size),
	    plot.margin=margin(t=30, r=0, b=10, l=5, unit="pt"))

    ### Calculate plot dimensions dependent on the sample count for spacing
    dpi <- 200
    width_px <- dpi * max(15, (9 + ((length(unique(meta.data$sample)) - 1) * 3)))
    height_px <- dpi * height

    ## Save the plot as PNG or SVG based on specified file format
    if (plot.file == "png") {
        png(file=file.name, width=width_px, height=height_px, units="px", res=dpi, type="cairo")
        try(print(p1_l1 + labs(title = title.plot)), silent=FALSE)
	dev.off()
    } else {
        p1_l1 <- p1_l1 + labs(title=title.plot)
        filename <- substr(file.name, 1, nchar(file.name) - 4)
	ggsave(
	    plot=p1_l1,
	    filename=paste0(filename, '.svg'),
	    width=length(unique(seurat.obj@meta.data[, split.by])) * 3,
	    height=height,
	    dpi=dpi
	)
    }
}

## Generates heatmap
heatmap_plot <- function(obj, dir, subset="FALSE", subset.size=50000,
			 assay.use="RNA", markers.csv="FALSE", metadata.vector) {
    cat("Value of dir: ", dir, "\n")
    cat("Value of subset: ", subset, "\n")
    cat("Value of subset.size: ", subset.size, "\n")
    cat("Value of assay.use: ", assay.use, "\n")
    cat("Value of markers.csv: ", markers.csv, "\n")

    ## Generate a list to hold the paths of the output plots
    heatmap_plots <- list()

    ## Generate the markers directory to hold the markers CSV files
    ## if they are getting generated in the initial run
    markers.dir <- paste0(dir, "/markers/")
    if (markers.csv == "FALSE") {
        dir.create(markers.dir)
	cat("Directory", markers.dir, "was generated.", "\n")
    }

    ## Subset the ingrated dataset if specified
    if (subset == "TRUE") {
	## Subsetting occurs when there are enough cells present in dataset
	print(nrow(obj@meta.data))
        if (nrow(obj@meta.data) > subset.size) {
            obj <- obj[, sample(colnames(obj), size=subset.size, replace=FALSE)]
	    print(nrow(obj@meta.data))
            cat("\nThe integrated dataset was subset.\n")
	} else {
            cat(paste0("\nThere were not enough cells present in the ",
		       "integrated dataset to subset.\n"))
	}
    } else {
        cat("\nThe integrated dataset was not subset.\n")
    }

    meta.data.list.plots.heatmap <- colnames(obj@meta.data)[
	                            colnames(obj@meta.data) %in% metadata.vector &
	                            rapply(obj@meta.data, function(x) length(unique(x))) < 16 &
			            rapply(obj@meta.data, function(x) length(unique(x))) > 1]
    if (assay.use == "SCT") {
        DefaultAssay(obj) <- "SCT"
    } else {
        DefaultAssay(obj) <- "RNA"
    }
    print(DefaultAssay(obj))
    obj <- NormalizeData(obj)

    for (meta in meta.data.list.plots.heatmap) {
        Idents(obj) <- meta
        csv.path <- paste0(markers.dir, meta, "_all_markers.csv")
	print(csv.path)
        if (markers.csv == "TRUE") {
	    cat("Marker CSV file from initial run located in" , markers.dir, "will be used.")
	    obj.markers <- read.csv2(csv.path, sep=",", row.names=1)
	} else {
	    if (assay.use == "SCT") {
                DefaultAssay(obj) <- "SCT"
	        obj <- PrepSCTFindMarkers(obj)
		obj.markers <- FindAllMarkers(obj, only.pos=TRUE, min.pct=0.05,
					      logfc.threshold=0.1, recorrect_umi=FALSE)
            } else {
                DefaultAssay(obj) <- "RNA"
	        obj.markers <- FindAllMarkers(obj, only.pos=TRUE, min.pct=0.05,
					      logfc.threshold=0.1)
	    }
            print(head(obj.markers, 10))
	    ## Write CSV file to markers directory
	    write.csv(obj.markers, file=csv.path)

	}

	if(nrow(unique(obj[[meta]])) < 5) {
	    marker.num <- 10
	} else if (nrow(unique(obj[[meta]]))) {
            marker.num <- 6
        } else {
            marker.num <- 4
	}
	print(marker.num)
        ## Identify the top markers within each cluster sorted by highest average log fold change
	top5 <- obj.markers %>% group_by(cluster) %>% top_n(n=marker.num, wt=avg_log2FC)
	feat.hm <- top5$gene
	print(feat.hm)

	obj <- ScaleData(obj, verbose=FALSE, assay=assay.use, features=feat.hm)

	## Generate heatmap plot
	dpi <- 200
	heatmap.path <- paste0(dir,"/cluster_heatmap_",meta,".png")
	plot_object <- DoHeatmap(obj, features=feat.hm, assay=assay.use, slot="scale.data", label=FALSE) + scale_fill_viridis()
	ggsave(filename=heatmap.path, plot=plot_object, dpi=200, width=20, height=12)

	## Add plot to group of files that will be moved into their own directory
	heatmap_plots <- append(heatmap_plots, heatmap.path)

    }
    return (heatmap_plots)
}

# Create cell type annotation plots for the report
create_dimplot <- function(obj, group) {
    dimplot <- DimPlot(obj, group.by=group, label=TRUE, pt.size=1.5,
		       raster=FALSE, repel=TRUE) +
               guides(color=guide_legend(ncol=1)) +
	       theme(legend.text=element_text(size=9)) +
	       theme(legend.key=element_rect(fill="snow2", color="white")) +
	       scale_color_discrete(labels=function(x) str_wrap(x, width=15))
    return (dimplot)
}

## Create spatial cell type annotation plots for the report
## TODO: The legend generated in SpatialDimPlot cannot be modified with
## guides() or scale_color_discrete() for unknown reasons; therefore,
## the legend from DimPlot is used as the legend for the SpatialDimPlot
## as a temporary solution.
create_spatial_dimplot <- function(obj, group, visium.type) {
    if (visium.type == "hd") {
        point.size.factor <- 10
    } else {
        point.size.factor <- 1.6
    }
    spatial.dimplot <- SpatialDimPlot(obj, group.by=group, label=FALSE, pt.size.factor=point.size.factor) +
	               theme(legend.position = "none") +
		       ggtitle(group) +
		       theme(plot.title=element_text(hjust=0.5, size=16))
    dimplot <- create_dimplot(obj, group)
    legend <- get_legend(dimplot)
    spatial.dimplot <- spatial.dimplot + plot_layout(guides="collect") + legend +
	               plot_layout(widths=c(0.82, 0.18))
    return (spatial.dimplot)
}
