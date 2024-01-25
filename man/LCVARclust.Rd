\name{LCVARclust}
\alias{LCVARclust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting Latent Class VAR Model
}
\description{
Function to fit a Latent Class VAR model with a given number of latent classes
}
\usage{
LCVARclust(Data, yVars, Time, ID, xContinuous = NULL,
           xFactor = NULL, Covariates = "equal-within-clusters",
           Clusters, Lags, smallestClN = 3,
           ICType = c("HQ", "SC"), RndSeed,
           Rand = 1, Rational = c(TRUE, FALSE),
           Initialization = NULL, SigmaIncrease = 10,
           it = 25, Conv = 1e-06, pbar = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
The data provided in a data.frame.
}
  \item{yVars}{
An integer vector specifying the position of the column(s) in dataframe \code{Data} that contain the endogenous variables (= the VAR time series)
}
  \item{Time}{
An integer specifying the position of the column in dataframe \code{Data} that contains the time point.
}
  \item{ID}{
An integer specifying the position of the column in dataframe \code{Data} that contains the ID variable for every participant.
}
  \item{xContinuous}{
An integer vector specifying the position of the column(s) in dataframe Data that contain the continuous exogenous variable(s), if present.
}
  \item{xFactor}{
An integer vector specifying the position of the column(s) in dataframe Data that contain the categorical exogenous variable(s), if present.
}
  \item{Covariates}{
Constraints on the parameters of the exogenous variable(s). So far only \code{Covariates = "equal-within-clusters"} can be specified.
}
  \item{Clusters}{
An integer vector specifying the numbers of latent classes (i.e., clusters) for which LCVAR models are to be calculated.
}
  \item{Lags}{
An integer specifying the number of VAR(p) lags to consider. Needs to be a subsequent sequence of subsequent integers.
}
  \item{smallestClN}{
An integer specifying the lowest number of individuals allowed in a cluster. When during estimation the crisp cluster membership of a cluster indicates less than \code{smallestClN} individuals, the covariance matrix and the posterior probabilities of cluster membership are reset.
}
  \item{ICType}{
The information criterion used to select the ideal model for a given number of clusters across all EM-starts and lag combinations. One of \code{c("HQ", "SC")}.
}
  \item{RndSeed}{Optional argument specifying the integer supplied to \code{set.seed()}, which guarantees reproducible results. If not specified, no seed is set.
}
  \item{Rand}{
The number of pseudo-random EM-starts used in fitting each possible model. For pseudo-random starts K individuals are randomly selected as cluster centres. Then individuals are partitioned into the cluster to which their ideographic VAR and ideographic covariate coefficients are closest.
}
  \item{Rational}{
Local, indicateing whether a rational EM-start should be used. Defaults to \code{Rational = TRUE}. Rational starts are based on the k-means partitioning of individuals’ ideographic VAR and ideographic covariate coefficients.
}
  \item{Initialization}{
An integer specifying the position of a column in dataframe \code{Data} that contains a guess at participants' cluster membership for a fixed number of clusters, if available.
}
  \item{SigmaIncrease}{
A numerical value specifying the value by which every element of Sigma will be increased when posterior probabilities of cluster memberships are reset.
}
  \item{it}{
An integer specifying the maximum number of EM-iterations allowed after every EM-start. After \code{it} EM-iterations an EM-start is forced to terminate.
}
  \item{Conv}{
A numerical value specifying the convergence criterion of the log likelihood to determine convergence of an EM-start. For details see Ernst et al. (2020) Inter-individual differences in multivariate time series: Latent class vector-autoregressive modelling.
}
  \item{pbar}{
If TRUE, a progress bar is shown. Defaults to \code{pbar = TRUE}.}
}
\details{
Some more details.
}

\value{
Describe everything in the output list ...
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{comp1}{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Ernst, A. F., Albers, C. J., Jeronimus, B. F., & Timmerman, M. E. (2020). Inter-individual differences in multivariate time-series. European Journal of Psychological Assessment.
}
\author{
Anja Ernst
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

\dontrun{

head(ExampleData)

LCVAR_out <- LCVARclust(Data = ExampleData, yVars = 1:4, Time = 6, ID = 5,
                        Covariates = "equal-within-clusters",
                        Clusters, Lags, smallestClN = 3,
                        ICType = "HQ", seme = 3, Rand = 2, Rational = TRUE,
                        SigmaIncrease = 10, it = 25, Conv = 1e-06, xContinuous = 7, xFactor = 8)


LCVAR_out$BestSolutionsPerCluster[[1]]$A


}


}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
