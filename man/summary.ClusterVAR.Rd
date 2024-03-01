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
An output object of the \code{ClusterVAR} function.
}
  \item{show}{
Indicate how models should be summarized. \code{show = "BPC"} compares models with different time lags for each number of clusters. For each number of clusters, the model with the best lag is selected and displayed in the output. The best lag is selected through the time-series information criterion selected with the argument \code{TS_criterion} (see below). \code{show = "GNC"} shows all models with different lags for a given number of clusters, this number of clusters is specified through \code{Number_of_Clusters} (see below). \code{show = "GNL"} shows for each number of clusters the model where all lags are fixed to a given number, this number of lags is specified through \code{Number_of_Lags} (see below). Out of these models, the best model in terms of the number of clusters is selected by the information criterion selected with the argument \code{global_criterion} (see below). Defaults to \code{show = "BPC"}.
}
\item{TS_criterion}{
The information criterion to select between different lags in a given cluster. The possible choices are \code{"SC"} and \code{"HQ"}. Defaults to \code{TS_criterion = "SC"}.
}
\item{global_criterion}{
The information criterion to select between different numbers of clusters for a given lag. The possible choices are \code{"BIC"} and \code{"ICL"}. Defaults to \code{global_criterion = "BIC"}.
}
\item{Number_of_Clusters}{
Specify the fixed number of clusters when specifying \code{show = "GNC"}. Defaults to \code{Number_of_Clusters = NULL}.
}
\item{Number_of_Lags}{
Specify the fixed number of lags when specifying \code{show = "GNL"}. Defaults to \code{Number_of_Lags = NULL}.
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




