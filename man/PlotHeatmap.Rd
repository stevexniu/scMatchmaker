% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{PlotHeatmap}
\alias{PlotHeatmap}
\title{Interaction Heatmap}
\usage{
PlotHeatmap(
  object,
  data = c("strength", "pvalue", "data"),
  which.cell = NULL,
  which.gene = NULL,
  aggregate = c("none", "cell", "gene"),
  counted = FALSE,
  selected = FALSE,
  scale = c("none", "cell", "gene"),
  scale_max = Inf,
  color = "RdBu",
  color.bins = 7,
  order = TRUE,
  order.dist = "minkowski",
  p = 1,
  filter = FALSE,
  low.cutoff = -Inf,
  high.cutoff = Inf,
  return.data = FALSE,
  Rowv = NA,
  Colv = NA,
  trace = "none",
  margins = c(10, 10),
  na.color = "black",
  key = TRUE,
  show.text = TRUE,
  ...
)
}
\arguments{
\item{object}{Matchmaker object.}

\item{data}{Data to plot. Default is interaction strength.
\itemize{
  \item strength, Interaction strength.
  \item pvalue, Interaction p value.
  \item data, expression data.
  }}

\item{which.cell}{Interacting cell type name (saperated by '|' i.e 'DC|Mac') or index to plot. Default is NULL.}

\item{which.gene}{Interacting gene name (saperated by '_' i.e 'CSF1_CSF1R') or index to plot. Default is NULL.}

\item{aggregate}{Type of aggregation used. Default is none.
\itemize{
  \item none, no aggregation.
  \item cell, aggregated by cell-cell interactions.
  \item gene, aggregated by ligand-receptor pairs.
  }}

\item{counted}{Aggregate interactions counts instead of strengths. Default is FALSE.}

\item{selected}{Use selected data if calculated. Default is FALSE.}

\item{scale}{Scale and center data. Default is none.
\itemize{
  \item none, no scaling.
  \item cell, scale data by cell (row).
  \item gene, scale data by gene (column).
  }}

\item{scale_max}{Maximum cutoff for scaled data. Default is Inf.}

\item{color}{Color palette name. Default is "RdBu", high values are in red and low values in blue. See details \code{\link[RColorBrewer]{brewer.pal}}.}

\item{color.bins}{Number of differen colors in pallete. Default is 7. See detail \code{\link[RColorBrewer]{brewer.pal}}.}

\item{order}{Whether to order the input data. Default is TRUE. See details \code{\link[seriation]{seriate}}.}

\item{order.dist}{Distance measure to be used. Default is "minkowski". See details \code{\link[stats]{dist}}.}

\item{p}{Power of Minkowski distance Default is 1, aka Manhattan distance. See details \code{\link[stats]{dist}}.}

\item{filter}{Whether to filter the input data based on low.cutoff and high.cutoff. Default is FALSE.}

\item{low.cutoff}{Lower cutoff bound below which will be removed. Default is negative infinity.}

\item{high.cutoff}{Upper cutoff bound above which will be removed. Default is positive infinity.}

\item{return.data}{Whether to return the data used to plot. Default is FALSE.}

\item{Rowv}{Whether to hierarchical cluster the rows. Default is NA. See details \code{\link[gplots]{heatmap.2}}.}

\item{Colv}{Whether to hierarchical cluster the columns Default is NA. See details \code{\link[gplots]{heatmap.2}}.}

\item{trace}{Whether to show the trace line. Defaults is none. See details \code{\link[gplots]{heatmap.2}}.}

\item{margins}{Plot margins. Defaults is c(10,10) See details \code{\link[gplots]{heatmap.2}}.}

\item{na.color}{Color for missing values (NAs). Defaults is black. See details \code{\link[gplots]{heatmap.2}}.}

\item{key}{Whether to show the key. Defaults is TRUE. See details \code{\link[gplots]{heatmap.2}}.}

\item{show.text}{Whether to show title, x-axis and y-axis text. Defaults is TRUE.}

\item{...}{Additioanl arguments passed to \code{\link[gplots]{heatmap.2}}.}
}
\value{
Plot a heatmap and optionally return the data used to generate the plot.
}
\description{
Plot a interaction heatmap.
}
\examples{
\dontrun{
 PlotHeatmap(object, data = "strength", aggregate = "row")
}
}
\concept{plot}
