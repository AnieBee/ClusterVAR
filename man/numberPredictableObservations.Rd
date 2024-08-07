\name{numberPredictableObservations}
\alias{numberPredictableObservations}

\title{
Determine the number of observations that can be predicted
}
\description{
\code{numberPredictableObservations} is a function to determine the number of observations in a given dataset that can be predicted based on the availability of previous observations, considering a specified time-lag.
}
\usage{
numberPredictableObservations(Data, yVars, Beep, Day = NULL, ID,
                              xContinuous = NULL, xFactor = NULL, Lags, ...)

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
Optional. An integer specifying the position of the column in dataframe \code{Data} that contains the variable that indicates the day of measurement. If \code{Day} is supplied here, measurements on the previous day are not used to predict measurements on the current day. Instead, the first \code{Lags} observations within each day are excluded from the calculation of VAR coefficients.
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
\item{Predictable observations per subject}{The number of predictable observations for each endogenous variable per subject, considering a specified time-lag.}
\item{Total predictable observations}{The total number of predictable observations summed over all subjects in the dataset for each endogenous variable, considering a specified time-lag.}
}
\author{
Anja Ernst & Jonas Haslbeck
}

\examples{
\donttest{
head(SyntheticData)

Obs <- numberPredictableObservations(Data = SyntheticData, yVars = 1:4,
                                      Beep = 9, Day = 10,  ID = 5, Lags = 1:3)

Obs

Obs$`Predictable observations per subject`$`1 Lag`
}

}


