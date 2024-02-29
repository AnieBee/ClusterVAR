\name{LCVAR}
\alias{LCVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting Latent Class VAR Model
}
\description{
Function to fit a Latent Class VAR model with a given number of latent classes
}
\usage{
LCVAR(Data, yVars, Beep, Day = NULL, ID,
           xContinuous = NULL, xFactor = NULL,
           Clusters, Lags, smallestClN = 3, RndSeed,
           Rand = 50, Rational = TRUE,
           Initialization = NULL, SigmaIncrease = 10,
           it = 50, Conv = 1e-05, pbar = TRUE,
           Covariates = "equal-within-clusters", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
The data provided in a data.frame.
}
  \item{yVars}{
An integer vector specifying the position of the column(s) in dataframe \code{Data} that contain the endogenous variables (= the VAR time series).
}
  \item{Beep}{
An integer specifying the position of the column in dataframe \code{Data} that contains the time point.
}
  \item{Day}{
Optional. An integer specifying the position of the column in dataframe \code{Data} that contains the variable that indicates on which day the measurement was made. If day is supplied here, the first measurment on each day will be removed from calculations for the VAR coefficients.
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
An integer or integer vector specifying the number of VAR(p) lags to consider. Needs to be a subsequent sequence of subsequent integers.
}
  \item{smallestClN}{
An integer specifying the lowest number of individuals allowed in a cluster. When during estimation the crisp cluster membership of a cluster indicates less than \code{smallestClN} individuals, the covariance matrix and the posterior probabilities of cluster membership are reset. Defaults to \code{smallestClN = 3}.
}
 \item{RndSeed}{
  Optional argument. An integer specifying the value supplied to \code{set.seed()}, which guarantees reproducible results. If not specified, no seed is set.
}
  \item{Rand}{
The number of pseudo-random EM-starts used in fitting each possible model. For pseudo-random starts K individuals are randomly selected as cluster centres. Then individuals are partitioned into the cluster to which their individual VAR and individual covariate coefficients are closest. Defaults to \code{Rand = 50}.
}
  \item{Rational}{
Logical, indicateing whether a rational EM-start should be used. Defaults to \code{Rational = TRUE}. Rational starts are based on the k-means partitioning of individuals’ ideographic VAR and ideographic covariate coefficients.
}
  \item{Initialization}{
Optional argument. An integer specifying the position of a column in dataframe \code{Data} that contains a guess at participants' cluster membership for a fixed number of clusters, if available.
}
  \item{SigmaIncrease}{
A numerical value specifying the value by which every element of Sigma will be increased when posterior probabilities of cluster memberships are reset.
}
  \item{it}{
An integer specifying the maximum number of EM-iterations allowed for every EM-start. After completing \code{it} EM-iterations, an EM-start is forced to terminate. Defaults to \code{it = 50}.
}
  \item{Conv}{
A numerical value specifying the convergence criterion of the log likelihood to determine convergence of an EM-start. For details see Ernst et al. (2020) Inter-individual differences in multivariate time series: Latent class vector-autoregressive modelling.  Defaults to \code{Conv = 1e-05}.
}
  \item{pbar}{
If TRUE, a progress bar is shown. Defaults to \code{pbar = TRUE}.}
  \item{Covariates}{
Constraints on the parameters of the exogenous variable(s). So far only \code{Covariates = "equal-within-clusters"} can be specified.}
  \item{...}{
Additional arguments passed to the function.}
}

\details{
Some more details.
}

\value{
An object of class 'ClusterVAR' providing several LCVAR models.
The details of the output components are as follows:
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{call}{The original function call.}
\item{Call}{A list of the arguments in the original function call.}
\item{Models}{All LCVAR models across all number of clusters, lag combinations, and number of starts. \code{Models[[a]][[b]][[c]]}  contains all information for the LCVAR model for the \code{a}th number of clusters that was specified in \code{Clusters}, for the \code{b}th combination of lag orders, on the \code{c}th EM-start. To find the ideal solution across all of them use \code{summary()} or \code{coef()}.}
\item{Runtime}{The runtime the function took to complete.}

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

head(Data)

LCVAR_out <- LCVAR(Data = Data,
                        yVars = 1:4, ID = 5, Beep = 9, Day = 10,
                        xContinuous = 7, xFactor = 8
                        Clusters = 1:2, Lags = 1:2, smallestClN = 3,
                        RndSeed = 3, Rand = 2, Rational = TRUE,
                        SigmaIncrease = 10, it = 25, Conv = 1e-05,
                        Covariates = "equal-within-clusters",)

summary(LCVAR_out)

summary(LCVARresult = LCVAR_out, show = "Given-a-number-of-clusters", Number_of_Clusters = 2)

summary(LCVARresult = LCVAR_out, show = "Best-overall")

coef(LCVAR_out, Lags = c(1))

coef(LCVAR_out, Lags = c(1, 1))

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