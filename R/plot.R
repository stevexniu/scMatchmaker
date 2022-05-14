#' @include utils.R
#' @importFrom graphics plot hist legend lines points text
#'
NULL

#' Interaction Spot Plot
#'
#' Generate a spot plot for cell-cell interactions similar to the bertin plot. 
#' The size or height of the dot or bar represents the interaction strength.
#' The color of the bars and dots represents the p values.
#'
#' @param object Matchmaker object.
#' @param which.cell Interacting cell type name (saperated by '|' i.e 'DC|Mac') or index to plot. Default is NULL.
#' @param which.gene Interacting gene name (saperated by '_' i.e 'CSF1_CSF1R') or index to plot. Default is NULL.
#' @param aggregate Type of aggregation used. Default is none.
#' \itemize{
#'   \item cell, aggregated by cell-cell interactions strengths (Default).
#'   \item gene, aggregated by ligand-receptor pairs strengths.
#'   \item none, no aggregation.
#'   }
#' @param counted Aggregate interactions counts instead of strengths. Default is FALSE.
#' @param selected Use selected data if calculated. Default is FALSE.
#' @param color Color palette name. Default is "RdBu" (reversed order), high values are in blue and low values in red See details \code{\link[RColorBrewer]{brewer.pal}}.
#' @param color.bins Number of differen colors in pallete. Default is 7. See detail \code{\link[RColorBrewer]{brewer.pal}}.
#' @param type Types of plot. Default is dotplot.
#' \itemize{
#'   \item dot, dotplot or bertin plot (Default).
#'   \item bar, barplot.
#'   }
#' @param order Whether to order the input data. Default is TRUE. See details \code{\link[seriation]{seriate}}.
#' @param order.dist Distance measure to be used. Default is "minkowski". See details \code{\link[stats]{dist}}.
#' @param p Power of Minkowski distance Default is 1, aka Manhattan distance. See details \code{\link[stats]{dist}}.
#' @param filter Whether to filter the input data based on low.cutoff and high.cutoff. Default is FALSE. 
#' @param low.cutoff Lower cutoff bound below which will be removed. Default is negative infinity.
#' @param high.cutoff Upper cutoff bound above which will be removed. Default is positive infinity.
#' @param dot.scale Scale the size of the dots. Default is 8.
#' @param plot.y.axis Whether to show y axis for each barplot. Default is TRUE.
#' @param fix.y.axis Whether to fix the all barplot ticks or use the relative ticks for each barplot. Default is TRUE. 
#' @param do.return Whether to return the data used to plot. Default is FALSE.
#' @return A ggplot object or a list with strength and p value data used to generate the plot. 
#' @export
#' @concept plot
#' @examples \dontrun{
#'  PlotSpot(object, which.cell = 1:10, which.gene = 1:10)
#' }
PlotSpot <- function(object, which.cell = NULL, which.gene = NULL, aggregate = c("none", "cell", "gene"), counted = FALSE, selected = FALSE, color = "RdBu", color.bins = 7, type = c("dot", "bar"),
                     order = FALSE, order.dist = "minkowski", p = 1, filter = FALSE, low.cutoff = -Inf, high.cutoff = Inf, dot.scale = 8, plot.y.axis = TRUE, fix.y.axis = TRUE, do.return = FALSE){
  aggregate <- match.arg(arg = aggregate)
  aggregate <- switch (aggregate,
                       cell = "row",
                       gene = "column",
                       none = "none")
  if(selected){
    if(length(x = object@selected) > 0){
      message("Using @selected data")
      spotPlot(interact.strength = object@selected$strength, interact.pval = object@selected$pvalue, which.cell = which.cell, which.gene = which.gene, aggregate = aggregate, counted = counted, color.bins = color.bins, color = color, 
               type = type, order = order, order.dist = order.dist, p = p, filter = filter, low.cutoff = low.cutoff, high.cutoff = high.cutoff, dot.scale = dot.scale, plot.y.axis = plot.y.axis, fix.y.axis = fix.y.axis, 
               do.return = do.return)
    } else {
      stop("@selected not found, please run Selecting(object, ...) first!", call. = FALSE)
    }
  } else {
    spotPlot(interact.strength = object@strength, interact.pval = object@pvalue, which.cell = which.cell, which.gene = which.gene, aggregate = aggregate, counted = counted, color.bins = color.bins, color = color, 
             type = type, order = order, order.dist = order.dist, p = p, filter = filter, low.cutoff = low.cutoff, high.cutoff = high.cutoff, dot.scale = dot.scale, plot.y.axis = plot.y.axis, fix.y.axis = fix.y.axis, 
             do.return = do.return)
  }
}

#' Interaction Heatmap
#'
#' Plot a interaction heatmap.
#'
#' @param object Matchmaker object.
#' @param data Data to plot. Default is interaction strength.
#' \itemize{
#'   \item strength, Interaction strength.
#'   \item pvalue, Interaction p value.
#'   \item data, expression data.
#'   }
#' @param which.cell Interacting cell type name (saperated by '|' i.e 'DC|Mac') or index to plot. Default is NULL.
#' @param which.gene Interacting gene name (saperated by '_' i.e 'CSF1_CSF1R') or index to plot. Default is NULL.
#' @param aggregate Type of aggregation used. Default is none.
#' \itemize{
#'   \item none, no aggregation.
#'   \item cell, aggregated by cell-cell interactions.
#'   \item gene, aggregated by ligand-receptor pairs.
#'   }
#' @param counted Aggregate interactions counts instead of strengths. Default is FALSE.
#' @param selected Use selected data if calculated. Default is FALSE.
#' @param scale Scale and center data. Default is none.
#' \itemize{
#'   \item none, no scaling.
#'   \item cell, scale data by cell (row).
#'   \item gene, scale data by gene (column).
#'   }
#' @param scale_max Maximum cutoff for scaled data. Default is Inf.
#' @param color Color palette name. Default is "RdBu", high values are in red and low values in blue. See details \code{\link[RColorBrewer]{brewer.pal}}.
#' @param color.bins Number of differen colors in pallete. Default is 7. See detail \code{\link[RColorBrewer]{brewer.pal}}.
#' @param order Whether to order the input data. Default is TRUE. See details \code{\link[seriation]{seriate}}.
#' @param order.dist Distance measure to be used. Default is "minkowski". See details \code{\link[stats]{dist}}.
#' @param p Power of Minkowski distance Default is 1, aka Manhattan distance. See details \code{\link[stats]{dist}}.
#' @param filter Whether to filter the input data based on low.cutoff and high.cutoff. Default is FALSE. 
#' @param low.cutoff Lower cutoff bound below which will be removed. Default is negative infinity.
#' @param high.cutoff Upper cutoff bound above which will be removed. Default is positive infinity.
#' @param return.data Whether to return the data used to plot. Default is FALSE.
#' @param Rowv Whether to hierarchical cluster the rows. Default is NA. See details \code{\link[gplots]{heatmap.2}}.
#' @param Colv Whether to hierarchical cluster the columns Default is NA. See details \code{\link[gplots]{heatmap.2}}.
#' @param trace Whether to show the trace line. Defaults is none. See details \code{\link[gplots]{heatmap.2}}.
#' @param margins Plot margins. Defaults is c(10,10) See details \code{\link[gplots]{heatmap.2}}.
#' @param na.color Color for missing values (NAs). Defaults is black. See details \code{\link[gplots]{heatmap.2}}.
#' @param key Whether to show the key. Defaults is TRUE. See details \code{\link[gplots]{heatmap.2}}.
#' @param show.text Whether to show title, x-axis and y-axis text. Defaults is TRUE.
#' @param ... Additioanl arguments passed to \code{\link[gplots]{heatmap.2}}.
#' @return Plot a heatmap and optionally return the data used to generate the plot.
#' @export
#' @concept plot
#' @examples \dontrun{
#'  PlotHeatmap(object, data = "strength", aggregate = "row")
#' }
PlotHeatmap <- function(object, data = c("strength","pvalue","data"), which.cell = NULL, which.gene = NULL, aggregate = c("none", "cell", "gene"), counted = FALSE, selected = FALSE, scale = c("none", "cell", "gene"), 
                        scale_max = Inf, color = "RdBu", color.bins = 7, order = TRUE, order.dist = "minkowski", p = 1, filter = FALSE, low.cutoff = -Inf, high.cutoff = Inf, 
                        return.data = FALSE, Rowv = NA, Colv = NA, trace = "none", margins = c(10, 10), na.color = "black", key = TRUE, show.text = TRUE, ...){
  data <- match.arg(arg = data)
  if(selected){
    if(length(x = object@selected) > 0){
      message("Using @selected data")
      data.use <- switch(data,
                         strength = object@selected$strength,
                         pvalues = object@selected$pvalue,
                         data = object@data)
    } else {
      stop("@selected not found, please run Selecting(object, ...) first!", call. = FALSE)
    }
  } else {
    data.use <- switch(data,
                       strength = object@strength,
                       pvalues = object@pvalue,
                       data = object@data)
  }
  scale <- match.arg(arg = scale)
  data.use <- switch (scale,
                      cell = Scaling(mat = data.use, scale_max = scale_max),
                      gene = t(x = Scaling(mat = t(x = data.use), scale_max = scale_max)),
                      none = as.matrix(x = data.use))
  aggregate <- match.arg(arg = aggregate)
  aggregate <- switch (aggregate,
                       cell = "row",
                       gene = "column",
                       none = "none")
  if(data == "data") aggregate = "none"
  heatmapPlot(data = data.use, which.cell = which.cell, which.gene = which.gene, aggregate = aggregate, counted = counted, color = color, color.bins = color.bins, order = order, order.dist = order.dist, p = p, filter = filter, 
              low.cutoff = low.cutoff, high.cutoff = high.cutoff, return.data = return.data, Rowv = Rowv, Colv = Colv, trace = trace, margins = margins, na.color = na.color, key = key, show.text = show.text, ...)
}

#' Interaction Histograms
#'
#' Plot a histogram for a pair of ligand and receptor.
#'
#' @param object Matchmaker object.
#' @param ligand Ligand name.
#' @param receptor Receptor name.
#' @param ligand.ident Cell identity to plot for ligand.
#' @param receptor.ident Cell identity to plot for receptor.
#' @param fix.breaks Whether to fix the \code{\link[graphics]{hist}} breaks using all the data. Default is TRUE.
#' @param nbins Number of bins to use for histogram. Default is 100.
#' @param freq Histogram of frequencies. Default is FALSE.
#' @param cols.use Two color code used for ligand and receptor. Default is pink and lightblue. 
#' @param alpha Color transparency. Default is 0.7.
#' @param ... Additional argument passing to \code{\link[graphics]{hist}}.	
#' @return Plot a histogram plot.
#' @export
#' @concept plot
#' @examples \dontrun{
#'  PlotHistogram(object, ligand = "CSF1", receptor = "CSF1R", 
#'  ligand.ident = "Tumor", receptor = "Macrophages")
#' }
PlotHistogram <- function(object, ligand, receptor, ligand.ident, receptor.ident, fix.breaks = TRUE, nbins = 100, freq = FALSE, cols.use = c("pink", "lightblue"), alpha = 0.7, ...){
  if(fix.breaks){
    breaks <- seq(from = min(object@data), to = max(object@data), length.out = nbins)
  } else{
    breaks <- NULL
  }
  histogramPlot(data = object@data, idents = object@annotation[,1], ligand = ligand, receptor = receptor, ligand.ident = ligand.ident, receptor.ident = receptor.ident, 
                breaks = breaks, nbins = nbins, freq = freq, cols.use = cols.use, alpha = alpha, ...)
}

#' Interaction Network Plot
#'
#' Generate a network plot of interactions.
#'
#' @param object Matchmaker object.
#' @param selected Use selected data if calculated. Default is FALSE.
#' @param which.cell Interacting cell type name (saperated by '|' i.e 'DC|Mac') or index to plot. Default is NULL.
#' @param which.gene Interacting gene name (saperated by '_' i.e 'CSF1_CSF1R') or index to plot. Default is NULL.
#' @param aggregate Type of aggregation used. Default is none.
#' \itemize{
#'   \item cell, aggregated by cell-cell interactions.
#'   \item gene, aggregated by ligand-receptor pairs.
#'   }
#' @param counted Aggregate interactions counts instead of strengths. Default is FALSE.
#' @param filter Whether to filter the input data based on low.cutoff and high.cutoff. Default is FALSE. 
#' @param low.cutoff Lower cutoff bound below which will be removed. Default is negative infinity.
#' @param high.cutoff Upper cutoff bound above which will be removed. Default is positive infinity.
#' @param edge.color Edge color. Default is black.
#' @param node.color Node color. Default is NULL.
#' @param node.size Node size. Default is 5
#' @param layout Layout used in the plot. Default is NULL.
#' @param directed Whether edges are directed. Default is TRUE.
#' @param bidirectional Whether edges are bidirectional. Default is TRUE.
#' @param arrows Whether to draw arrows. Default is TRUE.
#' @param louvain Whether to perform louvain clustering. Default is FALSE.
#' @param legend.posit Legend position. Default is topright.
#' @param legend.breaks Number of intervals to show. Default is 5.
#' @param legend.title Legend title. Default is 'Interaction Strength'.
#' @param ... Additioanl arguments passed to \code{\link[qgraph]{qgraph}}.
#' @return Return a network plot of interactions.
#' @export
#' @concept plot
#' @examples \dontrun{
#'  PlotNetwork(object, aggregate = "cell")
#' }
PlotNetwork <- function(object, selected = FALSE, which.cell = NULL, which.gene = NULL, aggregate = c("cell", "gene"), counted = FALSE, filter = FALSE, low.cutoff = -Inf, high.cutoff = Inf,
                        edge.color = "black", node.color = NULL, node.size = 5, layout = NULL, directed = TRUE, bidirectional = TRUE, arrows = TRUE, louvain = FALSE,
                        legend.posit = "topright", legend.breaks = 5, legend.title = "Interaction Strength", ...){
  if(selected){
    if(length(x = object@selected) > 0){
      message("Using @selected data")
      data <- object@selected$strength
    } else {
      stop("@selected not found, please run Selecting(object, ...) first!", call. = FALSE)
    }
  } else {
    data <- object@strength
  }
  aggregate <- match.arg(arg = aggregate)
  aggregate <- switch (aggregate,
                       cell = "row",
                       gene = "column")
  networkPlot(data = data, which.cell = which.cell, which.gene = which.gene, aggregate = aggregate, counted = counted, filter = filter, low.cutoff = low.cutoff, high.cutoff = high.cutoff, 
              edge.color = edge.color, node.color = node.color, node.size = node.size, layout = layout, directed = directed, bidirectional = bidirectional, arrows = arrows, louvain = louvain,
              legend.posit = legend.posit, legend.breaks = legend.breaks, legend.title = legend.title, ...)
}

#' Interaction Scatter Plot
#'
#' Generate a scatter plot of interactions.
#'
#' @param object Matchmaker object.
#' @param ident1 Cell identity 1 to plot.
#' @param ident2 Cell identity 2 to plot.
#' @param ligands Ligand name.
#' @param receptors Receptor name.
#' @param use_raw Logic, use raw data or ligand-receptor data. Default is TRUE.
#' @param add.lines Logic, draw lines between ligand-receptor pairs. Default is FALSE.
#' @param add.text Logic, add text for ligands and receptors. Default is TRUE.
#' @param background.col Background points color. Default is gray.
#' @param ligand.col Ligand points color. Default is red.
#' @param receptor.col Receptor points color. Default is blue.
#' @param ligand.pch Point shape for ligands. Default is 16 (solid circle).
#' @param receptor.pch Point shape for receptor. Default is 16 (solid circle).
#' @param point.cex Point size. Default is 1.
#' @param label.offset Text label offset. Default is 1. 
#' @param ligand.pos Ligand label position. Default is 2 (left). See \code{\link[graphics]{text}}.
#' @param receptor.pos Receptor label position. Default is 4 (right). See \code{\link[graphics]{text}}.
#' @param legend.pos Legend position. Default is 'topright'.
#' @param ... Additioanl arguments passed to \code{\link[base]{plot}}.
#' @return Return a scatter plot of interactions.
#' @export
#' @concept plot
#' @examples \dontrun{
#'  PlotScatter(object, ident1 = "Tumor", ident2 = "Immune", ligands = "CD274", receptors = "PDCD1")
#' }
PlotScatter <- function(object, ident1, ident2, ligands, receptors, use_raw = TRUE, add.lines = FALSE, add.text = TRUE, background.col = "grey", 
                        ligand.col = "red", receptor.col = "blue", ligand.pch = 16, receptor.pch = 16, point.cex = 1, label.offset = 1, 
                        ligand.pos = 2, receptor.pos = 4, legend.pos = "topright", ...){
  if(use_raw){
    data.use <- object@misc$raw_data
  } else {
    data.use <- object@data
  }
  scatterPlot(data = data.use, idents = object@annotation[,1], ident1 = ident1, ident2 = ident2, ligands = ligands, receptors = receptors, add.lines = add.lines, add.text = add.text, 
              background.col = background.col, ligand.col = ligand.col, receptor.col = receptor.col, ligand.pch = ligand.pch, receptor.pch = receptor.pch, point.cex = point.cex, 
              label.offset = label.offset, ligand.pos = ligand.pos, receptor.pos = receptor.pos, legend.pos = legend.pos, ...)
}

#' Generate Spot Plot
#'
#' Generate a spot barplot or dotplot.
#'
#' @keywords internal
#' @param interact.strength Interaction strength matrix.
#' @param interact.pval Interaction p-value matrix.
#' @param which.cell Interacting cell type name or index to plot. Default is NULL.
#' @param which.gene Interacting gene name or index to plot. Default is NULL.
#' @param aggregate Type of aggregation used.
#' \itemize{
#'   \item row, aggregated by row (Default).
#'   \item column, aggregated by column
#'   \item none, no aggregation.
#'   }
#' @param counted Aggregate interactions counts instead of strengths. Default is FALSE.
#' @param color Color palette name. Default is "RdBu" (reversed order), high values are in blue and low values in red See details \code{\link[RColorBrewer]{brewer.pal}}.
#' @param color.bins Number of differen colors in pallete. Default is 7. See detail \code{\link[RColorBrewer]{brewer.pal}}.
#' @param type Types of plot. Default is barplot.
#' \itemize{
#'   \item dot, dotplot or bertin plot (Default).
#'   \item bar, barplot.
#'   }
#' @param order Whether to order the input data. Default is TRUE. See details \code{\link[seriation]{seriate}}.
#' @param order.dist Distance measure to be used. Default is "minkowski". See details \code{\link[stats]{dist}}.
#' @param p Power of Minkowski distance Default is 1, aka Manhattan distance. See details \code{\link[stats]{dist}}.
#' @param filter Whether to filter the input data based on low.cutoff and high.cutoff. Default is FALSE. 
#' @param low.cutoff Lower cutoff bound below which will be removed. Default is negative infinity.
#' @param high.cutoff Upper cutoff bound above which will be removed. Default is positive infinity.
#' @param dot.scale Scale the size of the dots. Default is 8.
#' @param plot.y.axis Whether to show y axis for each barplot. Default is TRUE.
#' @param fix.y.axis Whether to fix the all barplot ticks or use the relative ticks for each barplot. Default is TRUE. 
#' @param do.return Whether to return the data used to plot. Default is FALSE.
#' @return A ggplot object or a list with strength and p value data used to generate the plot. 
#' @importFrom RColorBrewer brewer.pal
#' @importFrom reshape2 melt
#' @importFrom stats reformulate
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 ggplot geom_bar scale_fill_gradientn aes_string facet_grid xlab ylab theme theme_minimal element_blank element_text geom_point scale_size scale_color_gradientn 
spotPlot <- function(interact.strength, interact.pval, which.cell = NULL, which.gene = NULL, aggregate = c("none", "row", "column"), counted = FALSE, color = "RdBu", color.bins = 7, type = c("dot", "bar"),
                     order = FALSE, order.dist = "minkowski", p = 1, filter = FALSE, low.cutoff = -Inf, high.cutoff = Inf, dot.scale = 8, plot.y.axis = TRUE, fix.y.axis = TRUE, do.return = FALSE){
  if(!is.null(x = which.cell)){
    interact.strength <- as.matrix(x = interact.strength)[which.cell, ,drop = FALSE]
    interact.pval <-  as.matrix(x = interact.pval)[which.cell, ,drop = FALSE]
  } 
  if(!is.null(x = which.gene)){
    interact.strength <-  as.matrix(x = interact.strength)[,which.gene,drop = FALSE]
    interact.pval <- as.matrix(x = interact.pval)[,which.gene,drop = FALSE]
  } 
  aggregate <- match.arg(arg = aggregate)
  interact.strength <- switch(aggregate,
                              row = aggregateName(data = as.matrix(x = interact.strength), aggregate = "row", counted = counted),
                              column = aggregateName(data = as.matrix(x = interact.strength), aggregate = "column", counted = counted),
                              none = as.matrix(x = interact.strength))
  interact.pval <- as.matrix(x = interact.pval)
  pal.use <- brewer.pal(n = color.bins, name = color)
  if(filter){
    interact.strength <- filterMatrix(matrix = interact.strength, low.cutoff = low.cutoff, high.cutoff = high.cutoff)
    if(aggregate == "none") interact.pval <- interact.pval[rownames(x = interact.strength), colnames(x = interact.strength)]
  }
  if(nrow(x = interact.strength) == 1 | ncol(x = interact.strength) == 1) order <- FALSE
  if(order){
    mat.ordered <- orderMatrix(matrix = interact.strength, order.dist = order.dist, p = p, return.order = TRUE)
    interact.strength <- mat.ordered$ordered_matrix
    if(aggregate == "none") interact.pval <- permute(x = interact.pval, order = mat.ordered$order)
  } 
  if(aggregate != "none"){
    pal.use <- rev(x = pal.use)
    if(counted){
      strength.df <- data.frame(melt(data = interact.strength, varnames=c("interactions", "pairs"), value.name = "Counts"))
    } else {
      strength.df <- data.frame(melt(data = interact.strength, varnames=c("interactions", "pairs"), value.name = "Strength"))
    }
  } else {
    strength.df <- data.frame(melt(data = interact.strength, varnames=c("interactions", "pairs"), value.name = "Strength"))
    pval.df <- data.frame(melt(data = interact.pval, varnames=c("interactions", "pairs"), value.name = "pval"))
    strength.df$`p-values` <- pval.df$pval
  }
  type <- match.arg(arg = type)
  if(type == "bar"){
    if(aggregate != "none"){
      if(counted){
        p <- ggplot(data = strength.df, aes_string(x = 1, y = "Counts", fill = "Counts")) + geom_bar(stat = "identity") + 
          scale_fill_gradientn(colours = pal.use)
      } else {
        p <- ggplot(data = strength.df, aes_string(x = 1, y = "Strength", fill = "Strength")) + geom_bar(stat = "identity") + 
          scale_fill_gradientn(colours = pal.use)
      }
    } else {
      p <- ggplot(data = strength.df, aes_string(x = 1, y = "Strength", fill = "`p-values`")) + geom_bar(stat = "identity") + 
        scale_fill_gradientn(colours = pal.use)
    }
    if(fix.y.axis){
      p <- p + facet_grid(reformulate("interactions", "pairs"), switch = 'x') 
    } else {
      p <- p + facet_grid(reformulate("interactions", "pairs"), scales = "free_y",switch = 'x') + theme_cowplot()
    }
    if(aggregate == "row"){
      p <- p + xlab("Cell Types") + ylab("Cell Types")
    } else if(aggregate == "column"){
      p <- p + xlab("Interaction Partner A") + ylab("Interaction Partner B")
    } else {
      p <- p + xlab("Cell-Cell Interactions") + ylab("Interaction Strength")
    }
    p <- p + theme_minimal() + 
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            panel.grid = element_blank(),
            strip.text.y.right = element_text(angle = 0),
            strip.text.x = element_text(angle = 90),
            strip.background = element_blank())
    if(!plot.y.axis){
      p <- p + theme(axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.line.y = element_blank())
    }
  } else if(type == "dot"){
    strength.df$pairs <- factor(x = strength.df$pairs, levels = rev(x = levels(x = strength.df$pairs)), ordered = TRUE)
    p <- ggplot(data = strength.df, aes_string(x = "interactions", y = "pairs"))
    if(aggregate != "none"){
      if(counted){
        p <- p + geom_point(aes_string(size = "Counts", colour = "Counts"))
      } else {
        p <- p + geom_point(aes_string(size = "Strength", colour = "Strength"))
      }
    } else {
      p <- p + geom_point(aes_string(size = "Strength", colour = "`p-values`"))
    }
    p <- p + scale_size(range = c(0, dot.scale)) + scale_color_gradientn(colours = pal.use) + theme_cowplot() 
    if(aggregate == "row"){
      p <- p + xlab("Cell Types") + ylab("Cell Types")
    } else if(aggregate == "column"){
      p <- p + xlab("Interaction Partner A") + ylab("Interaction Partner B")
    } else {
      p <- p + xlab("Cell-Cell Interactions") + ylab("Ligand-Receptor Pairs")
    }
    p <- p + theme(axis.ticks.x=element_blank(), axis.text.x = element_text(angle = 90))    
  }
  if(do.return){
    return(list(ggplot = p, strength = interact.strength, pvalues = interact.pval))
  } else {
    return(p)
  }
}

#' Generate Heatmap
#'
#' Generate a heatmap plot.
#'
#' @keywords internal
#' @param data Data matrix to plot
#' @param which.cell Interacting cell type name or index to plot. Default is NULL.
#' @param which.gene Interacting gene name or index to plot. Default is NULL.
#' @param aggregate Type of aggregation used.
#' \itemize{
#'   \item row, aggregated by row (Default).
#'   \item column, aggregated by column
#'   \item none, no aggregation.
#'   }
#' @param counted Aggregate interactions counts instead of strengths. Default is FALSE.
#' @param color Color palette name. Default is "RdBu", high values are in red and low values in blue. See details \code{\link[RColorBrewer]{brewer.pal}}.
#' @param color.bins Number of differen colors in pallete. Default is 7. See detail \code{\link[RColorBrewer]{brewer.pal}}.
#' @param order Whether to order the input data. Default is TRUE. See details \code{\link[seriation]{seriate}}.
#' @param order.dist Distance measure to be used. Default is "minkowski". See details \code{\link[stats]{dist}}.
#' @param p Power of Minkowski distance Default is 1, aka Manhattan distance. See details \code{\link[stats]{dist}}.
#' @param filter Whether to filter the input data based on low.cutoff and high.cutoff. Default is FALSE. 
#' @param low.cutoff Lower cutoff bound below which will be removed. Default is negative infinity.
#' @param high.cutoff Upper cutoff bound above which will be removed. Default is positive infinity.
#' @param return.data Whether to return the data used to plot. Default is FALSE.
#' @param Rowv Whether to hierarchical cluster the rows. Default is NA. See details \code{\link[gplots]{heatmap.2}}.
#' @param Colv Whether to hierarchical cluster the columns Default is NA. See details \code{\link[gplots]{heatmap.2}}.
#' @param trace Whether to show the trace line. Defaults is none. See details \code{\link[gplots]{heatmap.2}}.
#' @param margins Plot margins. Defaults is c(10,10) See details \code{\link[gplots]{heatmap.2}}.
#' @param na.color Color for missing values (NAs). Defaults is black. See details \code{\link[gplots]{heatmap.2}}.
#' @param key Whether to show the key. Defaults is TRUE. See details \code{\link[gplots]{heatmap.2}}.
#' @param show.text Whether to show title, x-axis and y-axis text. Defaults is TRUE.
#' @param ... Additioanl arguments passed to \code{\link[gplots]{heatmap.2}}.
#' @return Plot a heatmap and optionally return the data used to generate the plot.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom gplots heatmap.2
heatmapPlot<- function(data, which.cell = NULL, which.gene = NULL, aggregate = c("none", "row", "column"), counted = FALSE, color = "RdBu", color.bins = 7, 
                       order = TRUE, order.dist = "minkowski", p = 1, filter = FALSE, low.cutoff = -Inf, high.cutoff = Inf, return.data = FALSE, 
                       Rowv = NA, Colv = NA, trace = "none", margins = c(10, 10), na.color = "black", key = TRUE, show.text = TRUE, ...){
  if(!is.null(x = which.cell)){
    if(length(x = which.cell) < 2) stop("which.cell cannot be single numbers.", call. = FALSE)
    data <- data[which.cell,]
  } 
  if(!is.null(x = which.gene)){
    if(length(x = which.gene) < 2) stop("which.gene cannot be single numbers.", call. = FALSE)
    data <- data[,which.gene]
  }
  pal.use <- rev(x = brewer.pal(n = color.bins, name = color))
  aggregate <- match.arg(arg = aggregate)
  aggregate.data <- switch(aggregate,
                           row = aggregateName(data = as.matrix(x = data), aggregate = "row", counted = counted),
                           column = aggregateName(data = as.matrix(x = data), aggregate = "column", counted = counted),
                           none = as.matrix(x = data))
  if(filter) aggregate.data <- filterMatrix(matrix = aggregate.data, low.cutoff = low.cutoff, high.cutoff = high.cutoff)
  if(is.null(x = dim(x = aggregate.data)) | any(dim(x = aggregate.data) < 2)) stop("Data has less than 2 rows and 2 columns after filtering", call. = FALSE)
  if(nrow(x = aggregate.data) == 1 | ncol(x = aggregate.data) == 1) order <- FALSE
  if(order) aggregate.data <- orderMatrix(matrix = aggregate.data, order.dist = order.dist, p = p, return.order = FALSE)
  if(all(aggregate.data == 0)) stop("Cannot plot when all entries are zero!", call. = FALSE)
  if(show.text){
    if(aggregate == "row"){
      xlab <- "Cell Type"
      ylab <- "Cell Type"
    } else if(aggregate == "column"){
      xlab <- "Interaction Partner A"
      ylab <- "Interaction Partner B"
    } else {
      xlab <- "Cell-Cell Interactions"
      ylab <- "Ligand-Receptor Pairs"
    }
    title <- ifelse(test = counted, yes = "Interaction Counts", no = "Interaction Strength")
    suppressWarnings(expr = heatmap.2(x = aggregate.data, Rowv = Rowv, Colv = Colv, trace = trace, col = pal.use, main = title, xlab = xlab, ylab = ylab, margins = margins, na.color = na.color, key = key, ...))
  } else {
    suppressWarnings(expr = heatmap.2(x = aggregate.data, Rowv = Rowv, Colv = Colv, trace = trace, col = pal.use, margins = margins, na.color = na.color, key = key, ...))
  }
  if(return.data) return(aggregate.data)
}

#' Generate Interaction Histograms
#'
#' Generate a histogram for a pair of ligand and receptor.
#'
#' @keywords internal
#' @param data Data matrix to plot
#' @param idents Cell type identity. 
#' @param ligand Ligand name.
#' @param receptor Receptor name.
#' @param ligand.ident Cell identity to plot for ligand.
#' @param receptor.ident Cell identity to plot for receptor.
#' @param breaks A vector giving the breakpoints between histogram cells. Default is NULL.
#' @param nbins Number of bins to use for histogram when breaks option is NULL. Default is 100.
#' @param freq Histogram of frequencies. Default is FALSE.
#' @param cols.use Two color code used for ligand and receptor. Default is pink and lightblue. 
#' @param alpha Color transparency. Default is 0.7
#' @param ... Additional argument passing to \code{\link[graphics]{hist}}.	
#' @return Plot a histogram.
#' @importFrom ggplot2 alpha
histogramPlot <- function(data, idents, ligand, receptor, ligand.ident, receptor.ident, breaks = NULL, nbins = 100, freq = FALSE, cols.use = c("pink", "lightblue"), alpha = 0.7, ...){
  if(length(x = ligand) > 1 | length(x = receptor) > 1 | length(x = ligand.ident) > 1 | length(x = receptor.ident) > 1) stop("Please provide only one ligand-receptor pair", call. = FALSE)
  ligand.data <- data[ligand, idents == ligand.ident]
  receptor.data <- data[receptor, idents == receptor.ident]
  if(is.null(x = breaks)) breaks <- seq(from = floor(x = min(c(ligand.data, receptor.data))), to = ceiling(x = max(c(ligand.data, receptor.data))), length.out = nbins)
  if(freq){
    counts_ligand <- hist(x = ligand.data, breaks = breaks ,plot = FALSE)$counts
    counts_receptor <- hist(x = receptor.data, breaks = breaks, plot = FALSE)$counts
    ylim.use <- max(c(counts_ligand, counts_receptor))
  } else {
    dens_ligand <- hist(x = ligand.data, breaks = breaks ,plot = FALSE)$density
    dens_receptor <- hist(x = receptor.data, breaks = breaks, plot = FALSE)$density
    ylim.use <- ceiling(x = max(c(dens_ligand, dens_receptor)))
  }
  hist(x = ligand.data, breaks = breaks, col = alpha(colour = cols.use[1], alpha = alpha), freq = freq, ylim = c(0, ylim.use), xlab = "Expression Level",
       main = paste(ligand,"-",receptor, " interaction in ", ligand.ident, "|", receptor.ident, sep = ""), border = alpha(colour = cols.use[1], alpha = alpha), ...)
  hist(x = receptor.data, breaks = breaks, add = TRUE, col = alpha(colour = cols.use[2], alpha = alpha), freq = freq, border = alpha(colour = cols.use[2], alpha = alpha))
  legend("topright", legend = c(paste(ligand, ligand.ident, sep = "-"), paste(receptor, receptor.ident, sep = "-")), fill = c("pink", "lightblue"), bty = "n")
}

#' Generate Network Plot
#'
#' Generate a network plot of interactions.
#'
#' @keywords internal
#' @param data Data matrix to plot
#' @param which.cell Interacting cell type name or index to plot. Default is NULL.
#' @param which.gene Interacting gene name or index to plot. Default is NULL.
#' @param aggregate Type of aggregation used.
#' \itemize{
#'   \item row, aggregated by row (Default).
#'   \item column, aggregated by column
#'   }
#' @param counted Aggregate interactions counts instead of strengths. Default is FALSE.
#' @param filter Whether to filter the input data based on low.cutoff and high.cutoff. Default is FALSE. 
#' @param low.cutoff Lower cutoff bound below which will be removed. Default is negative infinity.
#' @param high.cutoff Upper cutoff bound above which will be removed. Default is positive infinity.
#' @param edge.color Edge color. Default is black.
#' @param node.color Node color. Default is NULL.
#' @param node.size Node size. Default is 5
#' @param layout Layout used in the plot. Default is NULL.
#' @param directed Whether edges are directed. Default is TRUE.
#' @param bidirectional Whether edges are bidirectional. Default is TRUE.
#' @param arrows Whether to draw arrows. Default is TRUE.
#' @param louvain Whether to perform louvain clustering. Default is FALSE.
#' @param legend.posit Legend position. Default is topright
#' @param legend.breaks Number of intervals to show. Default is 5.
#' @param legend.title Legend title. Default is 'Interaction Strength'.
#' @param ... Additioanl arguments passed to \code{\link[qgraph]{qgraph}}.
#' @return Return a network plot of interactions.
#' @importFrom qgraph qgraph
#' @importFrom igraph as.igraph cluster_louvain
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
networkPlot <- function(data, which.cell = NULL, which.gene = NULL, aggregate = c("row", "column"), counted = FALSE, filter = FALSE, low.cutoff = -Inf, high.cutoff = Inf,
                        edge.color = "black", node.color = NULL, node.size = 5, layout = NULL, directed = TRUE, bidirectional = TRUE, arrows = TRUE, louvain = FALSE,
                        legend.posit = "topright", legend.breaks = 5, legend.title = "Interaction Strength", ...){
  if(!is.null(x = which.cell)) data <- data[which.cell, ,drop = FALSE]
  if(!is.null(x = which.gene)) data <- data[, which.gene,drop = FALSE]
  aggregate <- match.arg(arg = aggregate)
  aggregate.data <- switch(aggregate,
                           row = aggregateName(data = as.matrix(x = data), aggregate = "row", counted = counted),
                           column = aggregateName(data = as.matrix(x = data), aggregate = "column", counted = counted))
  if(filter){
    row.id <- apply(X = aggregate.data, MARGIN = 1, FUN = function(x) any(x > low.cutoff & x <= high.cutoff))
    col.id <- apply(X = aggregate.data, MARGIN = 2, FUN = function(x) any(x > low.cutoff & x <= high.cutoff))
    filter.id <- row.id | col.id
    aggregate.data <- aggregate.data[filter.id, filter.id, drop = FALSE]
  }
  if(nrow(x = aggregate.data) == 1) stop("Please select two different cell types.", call. = FALSE)
  if(louvain){
    g <- qgraph(input = aggregate.data, directed = FALSE, DoNotPlot = TRUE)
    g <- as.igraph(object = g, attributes = TRUE)
    clu <- cluster_louvain(graph = g)
  }
  if(is.null(x = node.color)){
    col.pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
    col.pals <- unique(x = unlist(x = mapply(FUN = brewer.pal, col.pals$maxcolors, rownames(x = col.pals))))
    color.idx <- 1:ncol(x = aggregate.data) %% length(x = col.pals)
    color.idx[which(x = color.idx == 0)] <- length(x = col.pals)
    if(louvain){
      node.color <- col.pals[clu$membership] 
    } else {
      node.color <- col.pals[color.idx]
    }
  } 
  g <- qgraph(input = aggregate.data, layout = layout, labels = colnames(x = aggregate.data), color = node.color, vsize = node.size, edge.color = edge.color, 
              color = node.color, directed = directed, bidirectional = bidirectional, arrows = arrows, ...)
  legend.use <- round(x = seq(from = min(aggregate.data), to = max(aggregate.data), length.out = legend.breaks), digits = 3)
  lwd <- round(x = seq(from = min(g$graphAttributes$Edges$width), to = max(g$graphAttributes$Edges$width), length.out = legend.breaks), digits = 3)
  col <- sort(x = g$graphAttributes$Edges$color, decreasing = TRUE)
  col <- col[round(x = seq(from = 1, to = length(x = col), length.out = legend.breaks), digits = 3)]
  legend(legend.posit, legend = legend.use, col = col, lwd = lwd, bty = "n", title = legend.title) 
}

#' Generate Scatter Plot
#'
#' Generate a scatter plot of interactions.
#'
#' @keywords internal
#' @param data Data matrix to plot
#' @param idents Cell type identity. 
#' @param ident1 Cell identity 1 to plot.
#' @param ident2 Cell identity 2 to plot.
#' @param ligands Ligand name.
#' @param receptors Receptor name.
#' @param add.lines Logic, draw lines between ligand-receptor pairs. Default is FALSE.
#' @param add.text Logic, add text for ligands and receptors. Default is TRUE.
#' @param background.col Background points color. Default is gray.
#' @param ligand.col Ligand points color. Default is red.
#' @param receptor.col Receptor points color. Default is blue.
#' @param ligand.pch Point shape for ligands. Default is 16 (solid circle).
#' @param receptor.pch Point shape for receptor. Default is 16 (solid circle).
#' @param point.cex Point size. Default is 1.
#' @param label.offset Text label offset. Default is 1. 
#' @param ligand.pos Ligand text position. Default is 2 (left). See \code{\link[graphics]{text}}.
#' @param receptor.pos Receptor text position. Default is 4 (right). See \code{\link[graphics]{text}}.
#' @param legend.pos Legend position. Default is 'topright'.
#' @param ... Additioanl arguments passed to \code{\link[base]{plot}}.
#' @return Return a scatter plot of interactions.
#' @importFrom Matrix rowMeans
scatterPlot <- function(data, idents, ident1, ident2, ligands, receptors, add.lines = FALSE, add.text = TRUE, background.col = "grey", 
                        ligand.col = "red", receptor.col = "blue", ligand.pch = 16, receptor.pch = 16, point.cex = 1, label.offset = 1, 
                        ligand.pos = 2, receptor.pos = 4, legend.pos = "topright", ...){
  cell1 <- which(idents == ident1)
  cell2 <- which(idents == ident2)
  cell1.data <- rowMeans(x = data[, cell1, drop = FALSE], na.rm = TRUE)
  cell2.data <- rowMeans(x = data[, cell2, drop = FALSE], na.rm = TRUE)
  ligand.points <- cbind(cell1.data[ligands], cell2.data[ligands])
  receptor.points<- cbind(cell1.data[receptors], cell2.data[receptors])
  plot(cell1.data, cell2.data, pch = 16, col = background.col, xlab = paste(ident1, "Expression"), ylab = paste(ident2, "Expression"), ...)
  points(ligand.points, pch = ligand.pch, col = ligand.col, cex = point.cex)
  points(receptor.points, pch = receptor.pch, col = receptor.col, cex = point.cex)
  if(add.text){
    text(ligand.points, pos = ligand.pos, labels = ligands, offset = label.offset, col= ligand.col)
    text(receptor.points, pos = receptor.pos, labels = receptors, offset = label.offset, col = receptor.col)
  }
  if(add.lines){
    for(i in 1:length(x = ligands)){
      lines(x = c(ligand.points[i, 1], receptor.points[i, 1]), y = c(ligand.points[i, 2], receptor.points[i, 2]), lty = 2, lwd = 2)
    }
  }
  legend(legend.pos, legend = c("Ligand", "Receptor"), pch = c(ligand.pch, receptor.pch), col = c(ligand.col, receptor.col))
}

# Circos not implemented 

# Sankey not implemented
