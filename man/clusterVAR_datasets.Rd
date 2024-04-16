\name{clusterVAR-datasets}
\title{Datasets included in the ClusterVAR package}

\alias{ExampleData}
\alias{SyntheticData}


\description{Internal clusterVAR datasets.}

\author{Anja Ernst & Jonas Haslbeck}

\details{

SyntheticData: Follows the structure where a 'Beep' variable and a 'Day' variable together indicate the time-point of the measurement. This data has missing observations that have been removed.

ExampleData: Follows the structure where a single 'Time' variable indicates the time-point of the measurement. This data is generated from a 2 cluster model wtih 1 lag for both clusters (i.e., a c(1, 1) Model).


}

\keyword{internal}

\examples{

\dontrun{

head(SyntheticData)

LCVAR_outExample1 <- LCVAR(Data = SyntheticData,
                        yVars = 1:4, ID = 5, Beep = 9, Day = 10,
                        xContinuous = 7, xFactor = 8,
                        Clusters = 1:2, Lags = 1, smallestClN = 3,
                        RndSeed = 3, Rand = 10, Rational = TRUE,
                        SigmaIncrease = 10, it = 25, Conv = 1e-05,
                        Covariates = "equal-within-clusters")
summary(LCVAR_outExample1)



head(ExampleData)
LCVAR_outExample2 <- LCVAR(Data = ExampleData,
                        yVars = 1:4, ID = 5, Beep = 6,
                        xContinuous = 7, xFactor = 8,
                        Clusters = 1:2, Lags = 1:2, smallestClN = 3,
                        RndSeed = 3, Rand = 2, Rational = TRUE,
                        SigmaIncrease = 10, it = 25, Conv = 1e-05,
                        Covariates = "equal-within-clusters")
summary(LCVAR_outExample2)
summary(object = LCVAR_outExample2, show = "GNC", Number_of_Clusters = 2)
coef(LCVAR_outExample2, Model = c(1, 1))



head(OneClusterDataSet)
LCVAR_outExample3 <- LCVAR(Data = OneClusterDataSet,
                        yVars = 1:4, ID = 5, Beep = 6,
                        xContinuous = 7, xFactor = 8,
                        Clusters = 1:2, Lags = 1, smallestClN = 3,
                        RndSeed = 3, Rand = 2, Rational = TRUE,
                        SigmaIncrease = 10, it = 25, Conv = 1e-05,
                        Covariates = "equal-within-clusters")
summary(LCVAR_outExample3)
coef(LCVAR_outExample3, Model = c(1))


}

}
