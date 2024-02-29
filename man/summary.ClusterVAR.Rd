\name{summary.ClusterVAR}
\alias{summary.ClusterVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of ClusterVAR objects
}
\description{
Takes the output of the ClusterVAR object and creates a small summary of the fitted model(s).
}
\usage{
\method{summary}{ClusterVAR}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An output object of the \code{ClusterVAR} function
}
  \item{show}{
Indicate which models to summarize. \code{show = "Best-per-number-of-clusters"} compares models with different numbers of clusters, where within cluster the mode with the best lag is selected. Within cluster, the best lag is selected with the information criterion selected with the argument \code{TS_criterion} (see below). \code{show = "Given-a-number-of-clusters"} compares models with different lags for a given number of clusters. And \code{show = "Best-overall"} shows the single model with the best lags (selected with \code{TS_criterion}) and number of clusters (selected with \code{global_criterion}). Defaults to \code{show = "Best-per-number-of-clusters"}.
}
\item{TS_criterion}{
The information criterion to select between different lags in a given cluster. The possible choices are \code{"SC"} and \code{"HQ"}. Defaults to \code{TS_criterion = "SC"}.
}
\item{global_criterion}{
The information criterion to select between different numbers of clusters, after selecting lags within each cluster using \code{TS_criterion}. The possible choices are \code{"BIC"} and \code{"ICL"}. Defaults to \code{global_criterion = "BIC"}.
}
\item{Number_of_Clusters}{
Specify the fixed number of clusters when specifying \code{show = "Given-a-number-of-clusters"}. Defaults to \code{Number_of_Clusters = NULL}.
}
  \item{\dots}{
Pass additional arguments.
}
}
\value{
\item{FunctionOutput}{Is a data frame containing summaries of the fitted models.}
}

\author{
Anja Ernst & Jonas Haslbeck
}

\examples{
# Add example
}




