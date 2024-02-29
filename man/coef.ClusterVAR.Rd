\name{coef.ClusterVAR}
\alias{coef.ClusterVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Show parameters of given Model
}
\description{
Extracts the parameters of a given model from the output object of the \code{LCVARclust} function.
}
\usage{
\method{coef}{ClusterVAR}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An output object of the \code{ClusterVAR} function
}
  \item{Lags}{
An integer vector specifying the mode for which parameters should be shown. For example, \code{Lags = c(1,2,1)} returns the parameters of a model with three clusters, in which the first cluster has one lag, the second cluster has 2 lags, and the third cluster has 1 lag.
}
  \item{\dots}{
Pass additional arguments.
}
}
\value{
STILL NEEDS TO BE DESCRIBED (can't do that now, because the function doesn't work for me)
}
\author{
Anja Ernst & Jonas Haslbeck
}
\examples{
# Add example later
}




