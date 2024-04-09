\name{plot.ClusterVAR}
\alias{plot.ClusterVAR}

\title{
Visualizing Model fit and information criteria
}
\description{
Creates a variety of plots summarizing fitted LCVAR models.
}
\usage{
plot.ClusterVAR(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An output object of the \code{ClusterVAR} function
}
  \item{show}{
Indicate summaries to plot. \code{show = "BPC"} compares models with different numbers of clusters, where within cluster the mode with the best lag is selected. Within cluster, the best lag is selected with the information criterion selected with the argument \code{TS_criterion} (see below). \code{show = "GNC"} compares models with different lags for a given number of clusters, specified with \code{Nclusters}. \code{show = "GNL"} compares models with different number of clusters, given the a number of lags specified with \code{Nlags}.

Alternatively, the VAR matrices of a specific model can be visualized. \code{show = specific} visualizes the VAR matrices across clusters for a given model which is specified with the argument \code{Model}. For example, \code{Model = c(1, 1, 1)} displays the VAR matrices for the model with three clusters, which each have a lag of 1. \code{show = specificDiff} displays the pairwise differences betwee clusters for a specific model, also specified with the argument\code{Model}. \code{show = specific} and \code{show = specificDef} are currently only implemented for models with lag equal 1.
}
\item{TS_criterion}{
The information criterion to select between different lags in a given cluster. The possible choices are \code{"SC"} and \code{"HQ"}. Defaults to \code{TS_criterion = "SC"}.
}
\item{global_criterion}{
The information criterion to select between different numbers of clusters, after selecting lags within each cluster using \code{TS_criterion}. The possible choices are \code{"BIC"} and \code{"ICL"}. Defaults to \code{global_criterion = "BIC"}.
}
\item{Number_of_Lags}{
Specify the fixed number of lags when specifying \code{show = "GNL"}. Defaults to \code{Number_of_Lags = NULL}.
\item{Number_of_Clusters}{
Specify the fixed number of clusters when specifying \code{show = "GNC"}. Defaults to \code{Number_of_Clusters = NULL}.
}
  \item{\dots}{
Pass additional arguments.
}
}

\details{
Creates different plots showing model fit or fit indices for a specified set of models.
}

\author{
Anja Ernst & Jonas Haslbeck
}

\examples{
# Add example
}
