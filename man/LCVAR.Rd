\name{LCVAR}
\alias{LCVAR}

\title{
Fitting Latent Class VAR Models
}
\description{
Function to fit a Latent Class VAR model with a given number of latent classes.
}
\usage{
LCVAR(Data, yVars, Beep, Day = NULL, ID,
           xContinuous = NULL, xFactor = NULL,
           Clusters, Lags, Center = FALSE,
           smallestClN = 3, Cores = 2,
           RndSeed = NULL, Rand = 50, Rational = TRUE,
           Initialization = NULL, SigmaIncrease = 10,
           it = 50, Conv = 1e-05, pbar = TRUE, verbose = TRUE,
           Covariates = "equal-within-clusters", ...)
}
\arguments{
  \item{Data}{
The data provided in a data.frame.
}
  \item{yVars}{
An integer vector specifying the position of the column(s) in dataframe \code{Data} that contain the endogenous variables (= the VAR time series).
}
  \item{Beep}{
An integer specifying the position of the column in dataframe \code{Data} that contains the time-point.
}
  \item{Day}{
Optional argument. An integer specifying the position of the column in dataframe \code{Data} that contains the variable that indicates the day of measurement. If \code{Day} is supplied here, measurements on the previous day are not used to predict measurements on the current day. Instead, the first \code{Lags} observations within each day are excluded from the calculation of VAR coefficients.
}
  \item{ID}{
An integer specifying the position of the column in dataframe \code{Data} that contains the ID variable for every participant.
}
  \item{xContinuous}{
Optional argument.  An integer vector specifying the position of the column(s) in dataframe Data that contain the continuous exogenous variable(s), if present. Exogenous variables are also known as covariates or as moderators for the within-person mean.
}
  \item{xFactor}{
Optional argument.  An integer vector specifying the position of the column(s) in dataframe Data that contain the categorical exogenous variable(s), if present. Exogenous variables are also known as covariates or as moderators for the within-person mean.
}
  \item{Clusters}{
An integer or integer vector specifying the numbers of latent classes (i.e., clusters) for which LCVAR models are to be calculated.
}
  \item{Lags}{
An integer or integer vector specifying the number of VAR(p) lags to consider. Needs to be a sequence of subsequent integers. The maximum number supported is \code{Lags = 3}.
}
  \item{Center}{
Logical, indicating whether the data (i.e., the endogenous variables) should be centered per person before calculations. If \code{Center = TRUE}, the differences in within-person means are removed from the data and the clustering is based only on similarity in VAR coefficients and (if exogenous variable(s) are specified) similarities in infleunces of exogenous variable(s). If \code{Center = FALSE}, the clustering is also based on similarities in within-person means. Defaults to \code{Center = FALSE}.}

\item{smallestClN}{
An integer specifying the lowest number of individuals allowed in a cluster. When during estimation the crisp cluster membership of a cluster indicates less than \code{smallestClN} individuals, the covariance matrix and the posterior probabilities of cluster membership are reset. Defaults to \code{smallestClN = 3}.
}
\item{Cores}{
A positive integer specifying the number of cores used to parallelize the computations. Specifying a high number of available cores can speed up computation. Defaults to \code{Cores = 2}.
}
 \item{RndSeed}{
  Optional argument. An integer specifying the value supplied to \code{set.seed()}, which guarantees reproducible results. If not specified, no seed is set.
}
  \item{Rand}{
The number of pseudo-random EM-starts used in fitting each possible model. For pseudo-random starts K individuals are randomly selected as cluster centres. Then individuals are partitioned into the cluster to which their individual VAR and individual covariate coefficients are closest. High numbers (e.g., 50 and above) ensure that a global optimum will be found, but will take longer to compute. Defaults to \code{Rand = 50}.
}
  \item{Rational}{
Logical, indicating whether a rational EM-start should be used in addition to the other EM-starts. Defaults to \code{Rational = TRUE}. Rational starts are based on the k-means partitioning of individuals’ ideographic VAR and ideographic covariate coefficients.
}
  \item{Initialization}{
Optional argument. An integer specifying the position of a column in dataframe \code{Data} that contains a guess at participants' cluster membership for a fixed number of clusters, if available. This initialization will be used as an additional EM-start.
}
  \item{SigmaIncrease}{
A numerical value specifying the value by which every element of Sigma will be increased when posterior probabilities of cluster memberships are reset. Defaults to \code{SigmaIncrease = 10}.
}
  \item{it}{
An integer specifying the maximum number of EM-iterations allowed for every EM-start. After completing \code{it} EM-iterations, an EM-start is forced to terminate. High numbers (e.g., 100 and above) ensure convergence, but will take longer to compute. Defaults to \code{it = 50}.
}
  \item{Conv}{
A numerical value specifying the convergence criterion of the log likelihood to determine convergence of an EM-start. For details see Ernst et al. (2020) Inter-individual differences in multivariate time series: Latent class vector-autoregressive modelling.  Defaults to \code{Conv = 1e-05}.
}
  \item{pbar}{
If \code{pbar = TRUE}, a progress bar is shown. Defaults to \code{pbar = TRUE}.}

\item{verbose}{If \code{verbose = FALSE}, output messages are limited. Additionally, the \code{pbar} argument is overridden, so the progress bar is not printed. Defaults to \code{verbose = TRUE}.}

  \item{Covariates}{
Constraints on the parameters of the exogenous variable(s). So far only \code{Covariates = "equal-within-clusters"} can be specified.}
  \item{...}{
Additional arguments passed to the function.}
}

\details{
This function estimates the latent class vector-autoregressive model to obtain latent classes (i.e., clusters) of individuals who are similar in VAR coefficients and (if specified) in within-person means and infleunces of exogenous variable(s).


\deqn{y_{i, t} =  w_{i, t} + \mu_{k} + \beta_{k} x_{i, t}}

\deqn{ w_{i, t} = (\sum_{a = 1}^{p} \Phi_{k, a} w_{i, t-a}) + u_{i, t}\qquad u_{i, t} \sim N(0, \Sigma_{k})}

Here \eqn{ \mu_{k}} represents an m x 1 vector that contains the cluster-wise conditional within-person mean for each y-variable in cluster k. \eqn{\beta_{k}} represents an m x q matrix that expresses the cluster-wise moderating influence of q exogenous variables (\eqn{x_{i, t}}) on the within-person means in cluster k.
\eqn{\Phi_{k, a}} represents an m×m matrix containing the cluster-wise VAR coefficients at lag a for cluster k. See the references below for details.
}

\value{
An object of class 'ClusterVAR' providing several LCVAR models.
The details of the output components are as follows:
\item{Call}{A list of arguments from the original function call.}
\item{All_Models}{All LCVAR models across all number of clusters, lag combinations, and number of EM-starts. \code{All_Models[[a]][[b]][[c]]}  contains all information for the LCVAR model for the \code{a}th number of clusters that was specified in \code{Clusters}, for the \code{b}th combination of lag orders, based on the combination of lags that was specified in \code{Lags}, on the \code{c}th EM-start. To find the ideal model across all of them use \code{summary()}, to view the coefficients of a given model, use \code{coef()}.}
\item{Runtime}{The runtime the function took to complete.}

}
\references{
Ernst, A. F., Albers, C. J., Jeronimus, B. F., & Timmerman, M. E. (2020). Inter-individual differences in multivariate time-series: Latent class vector-autoregressive modeling. European Journal of Psychological Assessment, 36(3), 482–491. \doi{10.1027/1015-5759/a000578}
}
\author{
Anja Ernst
}


\examples{

\donttest{

head(SyntheticData)
LCVAR_outExample1 <- LCVAR(Data = SyntheticData,
                           yVars = 1:4, ID = 5, Beep = 9, Day = 10,
                           xContinuous = 7, xFactor = 8,
                           Clusters = 1:2, Lags = 1,
                           Center = TRUE,
                           Cores = 2, # Adapt to local machine
                           RndSeed = 123, Rand = 1, it = 25)
summary(LCVAR_outExample1)
summary(object = LCVAR_outExample1, show = "GNL", Number_of_Lags = 1)
coef(LCVAR_outExample1, Model = c(1, 1))


head(ExampleData)
LCVAR_outExample2 <- LCVAR(Data = ExampleData,
                           yVars = 1:4, ID = 5, Beep = 6,
                           xContinuous = 7, xFactor = 8,
                           Clusters = 1:2, Lags = 1:2,
                           Center = FALSE,
                           Cores = 2, RndSeed = 123,
                           Rand = 1,
                           it = 25, Conv = 1e-05)
summary(LCVAR_outExample2)
summary(object = LCVAR_outExample2, show = "GNL", Number_of_Lags = 1)
summary(object = LCVAR_outExample2, show = "GNC", Number_of_Clusters = 2)
coef(LCVAR_outExample2, Model = c(1, 1))
plot(LCVAR_outExample2, show = "specific", Model = c(1, 1))


LCVAR_outExample3 <- LCVAR(Data = ExampleData,
                           yVars = c("Item1", "Item2","Item3", "Item4"),
                           ID = "Person", Beep = "Timepoint",
                           xContinuous = "ContiniousVariable",
                           xFactor = "CategoricalVariable",
                           Clusters = 1:2, Lags = 1:2,
                           Center = FALSE,
                           Cores = 2, RndSeed = 123,
                           Rand = 1,
                           it = 25, Conv = 1e-05)
plot(LCVAR_outExample3, show = "GNL", Number_of_Lags = 1)
plot(LCVAR_outExample3, show = "GNC", Number_of_Clusters = 2)
}

}

