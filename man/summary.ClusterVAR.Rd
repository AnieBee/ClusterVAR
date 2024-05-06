\name{summary.ClusterVAR}
\alias{summary.ClusterVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of ClusterVAR objects
}
\description{
Takes the output of the LCVAR function and creates a small summary of the fitted model(s).
}
\usage{
\method{summary}{ClusterVAR}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An output object of the \code{LCVAR} function.
}
  \item{show}{
Indicate how models should be summarized, the possible choices are \code{"BPC"}, \code{"GNC"}and \code{"GNL"}. \code{show = "BPC"} compares models with different time lags for each number of clusters. If \code{show = "BPC"}, for each number of clusters the model with the best lag is selected and displayed in the output. The best lag is selected through the time-series information criterion specified with the argument \code{TS_criterion} (see below). \code{show = "GNC"} shows all models with different lags for a given number of clusters, this number of clusters is specified through \code{Number_of_Clusters} (see below). \code{show = "GNL"} shows for each number of clusters the model where all lags are fixed to a given number, this number of lags is specified through \code{Number_of_Lags} (see below). Out of these models, the best model in terms of the number of clusters is selected by the information criterion selected with the argument \code{global_criterion} (see below). Defaults to \code{show = "BPC"}.
}
\item{TS_criterion}{
The information criterion to select the best model between models with a different number of lags but with the same number of clusters. The possible choices are \code{"SC"} and \code{"HQ"}. Defaults to \code{TS_criterion = "SC"}.
}
\item{global_criterion}{
The information criterion to select the best model between models with different numbers of clusters but with the same number of lags. The possible choices are \code{"BIC"} and \code{"ICL"}. Defaults to \code{global_criterion = "BIC"}.
}
\item{Number_of_Clusters}{
An integer. Specify the fixed number of clusters when using \code{show = "GNC"}. Defaults to \code{Number_of_Clusters = NULL}.
}
\item{Number_of_Lags}{
An integer. Specify the fixed number of lags when using \code{show = "GNL"}. Defaults to \code{Number_of_Lags = } the lowest number of lags specified in object.
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

\references{
Hamilton, J. (1994), Time Series Analysis, Princeton University Press, Princeton.

Hannan, E. J. and B. G. Quinn (1979), The determination of the order of an autoregression, Journal of the Royal Statistical Society.

Lütkepohl, H. (2006), New Introduction to Multiple Time Series Analysis, Springer, New York.

Quinn, B. (1980), Order determination for a multivariate autoregression, Journal of the Royal Statistical Society.

Biernacki, C., Celeux, G., & Govaert, G. (2000). Assessing a mixture model for clustering with
the integrated completed likelihood. IEEE Transactions on Pattern Analysis and Machine
Intelligence.

Schwarz, G. (1978). Estimating the Dimension of a Model. The Annals of Statistics.


}


\seealso{\code{plot.ClusterVAR()}, \code{coef.ClusterVAR()}}



\examples{

\dontrun{
LCVAR_outExample <- LCVAR(Data = ExampleData,
                          yVars = 1:4, ID = 5, Beep = 6,
                          xContinuous = 7, xFactor = 8,
                          Clusters = 1:2, Lags = 1:2,
                          Cores = 2, RndSeed = 3,
                          Rand = 2, it = 25)

summary(LCVAR_outExample)
summary(object = LCVAR_outExample, show = "GNL", Number_of_Lags = 1)
summary(object = LCVAR_outExample, show = "GNL", Number_of_Lags = 1, global_criterion = "ICL")
summary(object = LCVAR_outExample, show = "GNC", Number_of_Clusters = 2)
summary(object = LCVAR_outExample, show = "GNC", Number_of_Clusters = 2, TS_criterion = "HQ")

}
}



