\name{plot.ClusterVAR}
\alias{plot.ClusterVAR}

\title{
Visualizing Model fit and information criteria
}
\description{
Creates a variety of plots summarizing fitted LCVAR models.
}
\usage{
\method{plot}{ClusterVAR}(x,
    show,
    Number_of_Clusters = NULL,
    Number_of_Lags = NULL,
    Model = NULL,
    ...)
}

\arguments{
  \item{x}{
An output object of the \code{ClusterVAR} function
}
  \item{show}{
Indicate summaries to plot. \code{show = "GNC"} compares models with different lags for a given number of clusters, specified with \code{Number_of_Clusters}. \code{show = "GNL"} compares models with different number of clusters, given the a number of lags specified with \code{Number_of_Lags}.

Alternatively, the VAR matrices of a specific model can be visualized. \code{show = specific} visualizes the VAR matrices across clusters for a given model which is specified with the argument \code{Model}. \code{show = specificDiff} displays the pairwise differences between clusters for a specific model, also specified with the argument \code{Model}. \code{show = specific} and \code{show = specificDiff} are currently only implemented for models with 1 lag.
}

\item{Number_of_Clusters}{
An integer. Specify the fixed number of clusters when using \code{show = "GNC"}. Defaults to \code{Number_of_Clusters = NULL}
}
\item{Number_of_Lags}{
An integer. Specify the fixed number of lags when using \code{show = "GNL"}. Defaults to \code{Number_of_Lags =} the lowest number of lags specified in object.
}
\item{Model}{
An integer vector. Specify when using \code{show = "specific"} or \code{show = specificDef}. Indicates the model for which coefficients should be plotted. For example, \code{Model = c(1,1,1)} plots a model with three clusters and each cluster has 1 lag. Defaults to \code{Model = NULL}.
}
  \item{\dots}{
Pass additional arguments.
}
}

\details{
Creates different plots showing either a fitted LCVAR model or fit indices for a specified set of LCVAR models.
}

\author{
Anja Ernst & Jonas Haslbeck
}

\examples{
\donttest{
LCVAR_outExample <- LCVAR(Data = ExampleData,
                          yVars = 1:4, ID = 5, Beep = 6,
                          xContinuous = 7, xFactor = 8,
                          Clusters = 1:2, Lags = 1:2,
                          Center = FALSE, Cores = 2,
                          RndSeed = 3, Rand = 2,
                          it = 25)

plot(LCVAR_outExample, show = "GNL", Number_of_Lags = 1)
plot(LCVAR_outExample, show = "GNC", Number_of_Clusters = 2)
plot(LCVAR_outExample, show = "specific", Model = c(1, 1))
plot(LCVAR_outExample, show = "specific", Model = c(1, 1), labels = c("A", "B", "C","D"))
plot(LCVAR_outExample, show = "specificDiff", Model = c(1, 1))
}
}

