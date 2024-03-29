% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{spotPlot}
\alias{spotPlot}
\title{Generate Spot Plot}
\usage{
spotPlot(
  interact.strength,
  interact.pval,
  which.cell = NULL,
  which.gene = NULL,
  aggregate = c("none", "row", "column"),
  counted = FALSE,
  color = "RdBu",
  color.bins = 7,
  type = c("dot", "bar"),
  order = FALSE,
  order.dist = "minkowski",
  p = 1,
  filter = FALSE,
  low.cutoff = -Inf,
  high.cutoff = Inf,
  dot.scale = 8,
  plot.y.axis = TRUE,
  fix.y.axis = TRUE,
  do.return = FALSE
)
}
\arguments{
\item{interact.strength}{Interaction strength matrix.}

\item{interact.pval}{Interaction p-value matrix.}

\item{which.cell}{Interacting cell type name or index to plot. Default is NULL.}

\item{which.gene}{Interacting gene name or index to plot. Default is NULL.}

\item{aggregate}{Type of aggregation used.
\itemize{
  \item row, aggregated by row (Default).
  \item column, aggregated by column
  \item none, no aggregation.
  }}

\item{counted}{Aggregate interactions counts instead of strengths. Default is FALSE.}

\item{color}{Color palette name. Default is "RdBu" (reversed order), high values are in blue and low values in red See details \code{\link[RColorBrewer]{brewer.pal}}.}

\item{color.bins}{Number of differen colors in pallete. Default is 7. See detail \code{\link[RColorBrewer]{brewer.pal}}.}

\item{type}{Types of plot. Default is barplot.
\itemize{
  \item dot, dotplot or bertin plot (Default).
  \item bar, barplot.
  }}

\item{order}{Whether to order the input data. Default is TRUE. See details \code{\link[seriation]{seriate}}.}

\item{order.dist}{Distance measure to be used. Default is "minkowski". See details \code{\link[stats]{dist}}.}

\item{p}{Power of Minkowski distance Default is 1, aka Manhattan distance. See details \code{\link[stats]{dist}}.}

\item{filter}{Whether to filter the input data based on low.cutoff and high.cutoff. Default is FALSE.}

\item{low.cutoff}{Lower cutoff bound below which will be removed. Default is negative infinity.}

\item{high.cutoff}{Upper cutoff bound above which will be removed. Default is positive infinity.}

\item{dot.scale}{Scale the size of the dots. Default is 8.}

\item{plot.y.axis}{Whether to show y axis for each barplot. Default is TRUE.}

\item{fix.y.axis}{Whether to fix the all barplot ticks or use the relative ticks for each barplot. Default is TRUE.}

\item{do.return}{Whether to return the data used to plot. Default is FALSE.}
}
\value{
A ggplot object or a list with strength and p value data used to generate the plot.
}
\description{
Generate a spot barplot or dotplot.
}
\keyword{internal}
