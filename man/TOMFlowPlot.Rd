% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overlaps.R
\name{TOMFlowPlot}
\alias{TOMFlowPlot}
\title{TOMFlowPlot}
\usage{
TOMFlowPlot(
  WGCNAlist,
  networks,
  toms,
  genes_to_label,
  alpha = 0.1,
  color = "black",
  width = 0.05
)
}
\arguments{
\item{WGCNAlist}{list of WGCNA objects}

\item{networks}{list of network names of length 2}

\item{toms}{a list of TOM distance objects of length 2}

\item{genes_to_label}{genes to label across two networks}

\item{alpha}{alpha of flows}

\item{color}{color of flows}

\item{width}{width of the strata}
}
\value{
a ggplot object
}
\description{
Plots a sankey flow diagram showing the movement of genes from one WGCNA
to another WGCNA. Uses the ggalluvial framework.
}
\author{
Dario Tommasini
}
