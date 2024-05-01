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
An output object of the \code{LCVAR} function.
}
  \item{Model}{
An integer vector specifying the model for which coefficients should be shown. For example, \code{Model = c(1,1,2)} returns the coefficients of a model with three clusters (i.e., latent classes), in which the first cluster has one lag, the second cluster has 1 lag, and the third cluster has 2 lags. The vector length indicates the number of clusters in the model of interest.
}
  \item{\dots}{
Pass additional arguments.
}
}
\value{

\item{Lags}{An integer vector specifying for which model the coefficients are shown. For example, \code{Lags = c(1,1,2)} indicates a model with three clusters, in which the first cluster has one lag, the second cluster has 1 lag, and the third cluster has 2 lags. }

\item{Classification}{The crisp classification for each individual into a cluster. Crisp classifications are made based on an individual's  modal cluster membership probabilities.}

\item{VAR_coefficients}{The cluster-wise vector-autoregressive coefficients for each cluster.}

\item{Exogenous_coefficients}{The cluster-wise exogenous coefficients for each cluster. The first column in each array indicates the (conditional) within-person mean in that cluster. If exogenous variable(s) were specified, the other columns indicate the influences of the exogenous variable(s) in that cluster.}

\item{Sigma}{The cluster-wise innovation covariance matrix for each cluster.}

\item{Proportions}{The mixing proportions for each cluster. These can be considered as the proportion of individuals belonging to the respective cluster.}

\item{PredictableTimepoints}{The total number of time-points in this dataset that could be predicted because the previous time-point(s) were observed. In case of unequal lags across different clusters, the number of time-points for each person are weighted by their posterior cluster-membership probability.  See also \code{numberPredictableObservations()}.}

\item{Converged}{A logical value that indicates whether this model converged.}
}

\author{
Anja Ernst & Jonas Haslbeck
}
\examples{
LCVAR_outExample <- LCVAR(Data = ExampleData,
                           yVars = 1:4, ID = 5, Beep = 6,
                           xContinuous = 7, xFactor = 8,
                           Clusters = 2, Lags = 1:2, smallestClN = 3,
                           Cores = 1, RndSeed = 3, Rand = 2,
                           it = 25)

coef(LCVAR_outExample, Model = c(1, 1))
coef(LCVAR_outExample, Model = c(2, 2))
}




