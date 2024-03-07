\name{coef.ClusterVAR}
\alias{coef.ClusterVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Show coefficients of given Model
}
\description{
Extracts the coefficients of a given model from the output object of the \code{LCVAR} function.
}
\usage{
\method{coef}{ClusterVAR}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An output object of the \code{LCVAR} function
}
  \item{Model}{
An integer vector specifying the model for which coefficients should be shown. For example, \code{Lags = c(1,1,2)} returns the coefficients of a model with three clusters, in which the first cluster has one lag, the second cluster has 1 lag, and the third cluster has 2 lags.  The vector length indicates the number of clusters in the model of interest.
}
  \item{\dots}{
Pass additional arguments.
}
}
\value{
STILL NEEDS TO BE DESCRIBED (can't do that now, because the function doesn't work for me)

\item{PredictableTimepoints}{The total number of time-points in this dataset that could be predicted because the previous time-point(s) were observed. In case of unequal lags across different clusters, the number of time-points for each person are weighted by posterior cluster-membership probability. }
}
\author{
Anja Ernst & Jonas Haslbeck
}
\examples{
# Add example later
}




