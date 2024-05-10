\name{numberPredictableObservations}
\alias{numberPredictableObservations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determine the number of observations that can be predicted
}
\description{
\code{numberPredictableObservations} is a function to determine the number of observations in a given dataset that can be predicted based on previous observations, considering a specified time-lag.
}
\usage{
numberPredictableObservations(Data, yVars, Beep, Day = NULL, ID, Lags, ...)

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
An integer specifying the position of the column in dataframe \code{Data} that contains the time-point.
}
  \item{Day}{
Optional. An integer specifying the position of the column in dataframe \code{Data} that contains the variable that indicates the day of measurement. If \code{Day} is supplied here, measurements on the previous day are not used to predict measurements on the current day. Instead, the first \code{Lags} observations within each day are excluded from the calculation of VAR coefficients.
}
  \item{ID}{
An integer specifying the position of the column in dataframe \code{Data} that contains the ID variable for every participant.
}
  \item{Lags}{
An integer or integer vector specifying the number of VAR(p) lags to consider. Needs to be a sequence of subsequent integers. The maximum number supported is \code{Lags = 3}.
}
  \item{...}{
Additional arguments passed to the function.}
}

\details{
This function determines the number of observations in a given dataset that can be predicted based on previous observations. For instance, in a lag-1 model, if an observation is missing, the observation at the next time-point cannot be predicted. Similarly, in a lag-2 model, if an observation is missing, the observations at the next two time-points cannot be predicted. The output gives the number of predictable observations for each of the endogenous variables that was specified under \code{yVars}. The number of predictable observations is the same for all endogenous variables.
}

\value{
The total number of predictable observations for each endogenous variable, considering a specified time-lag.
}
\author{
Anja Ernst & Jonas Haslbeck
}

\examples{

\dontrun{

head(SyntheticData)

Obs <- numberPredictableObservations(Data = SyntheticData, yVars = 1:4,
                                      Beep = 9, Day = 10,  ID = 5, Lags = 1:3)

Obs
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
