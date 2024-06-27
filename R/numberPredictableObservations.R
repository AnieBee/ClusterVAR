

numberPredictableObservations <- function(Data,
                                          yVars,
                                          Beep,
                                          Day = NULL,
                                          ID,
                                          Lags,
                                          ...)
    # If each measurement is done on the same day, don't specify day but only specify beep.
    # Only if the first measurment on each day should be removed from calculations, the
    # Data = data frame containing one or several columns for: yVars (position of columns containing endogenous VAR process variables in dataframe Data),
{
    # ------ Computing Some Aux Vars ------
    LowestLag = min(Lags)
    HighestLag = max(Lags)


    # ------ Input Checks ------
    # Cluster Search Sequence

    # Lag Search Sequence
    if (missing(Lags))
        stop("Specify the sequence of number of lags to search.")
    if (!all(Lags == round(Lags)))
        stop("Lags need to be specified as integers.")
    if (!all((Lags[-1] - Lags[-length(Lags)] == 1)))
        stop("Lags need to be specified as subsequent integers.")

    # Checks
    stopifnot(HighestLag <= 3) # highest lag number that is allowed is 3
    stopifnot(LowestLag & HighestLag)
    stopifnot(HighestLag >= LowestLag)
    stopifnot(length(yVars) > 1) # So far only multivariate time-series are implemented
    stopifnot(is.numeric(Data[, Beep]))
    if (!is.null(Day))
        stopifnot(is.numeric(Data[, Day]))

    #call <- match.call()

    PreviousSol = TRUE # Use solution of previous Lags as a start

    # Y has to be numeric, X can be numeric or factor
    # Y is data of form m \times (sum^N (nObs) )
    # ID = has to be factor


    # Remove rows with NA values
    Data <- Data[stats::complete.cases(Data[, c(yVars, Beep, Day)]), ] # only consider relevant columns to avoid issues with irrelevant NA-structures in other

    ##### Preprocessing of Data Set #####--------------------
    Data = as.data.frame(Data)
    if (is.null(Day)) {
        Data = Data[order(Data[, ID], Data[, Beep]),]
    } else{
        Data = Data[order(Data[, ID], Data[, Day], Data[, Beep]),]
    }
    # order Data according to ID, make sure an individual's observations occur one after another with first obs first, second second etc
    # observations have to occur ascending in time



    # Endogenous Variables #-------------------
    Y = t(as.matrix(Data[, yVars]))
    nDepVar = dim(Y)[1]



    # ID indicator Variable & TS length for every person #------------------------------
    Pers = as.character(Data[, ID])
    pers = unique(Pers)
    N = length(unique(Pers)) # length(levels(Pers))

    nObs = rep(0, N)
    for (i in 1:N)
    {
        # determine the total number of Observations (nObs) per pers
        nObs[i] = length(which(Pers == pers[i]))
    }

    stopifnot(identical(rep(pers, c(nObs)), as.character(Data[, ID]))) # should be superflous is a check nObs is correct
    stopifnot(is.numeric(Y))


    # Whole Sample of Obs
    PersStart = cumsum(c(0, nObs[-length(nObs)])) + 1 # Start of individual time series for every pers in Y or W
    PersEnd = cumsum(nObs) # end of individual time series of every pers, last obs of every pers in Y or W


    ### ------------ whether an observation in Y can be predicted (PredictableObs) ------------------
    PredictableObs = vector("list", HighestLag)  # from 1:HighestLag but only LowestLag:HighestLag elements are filled
    for (lagRunner in LowestLag:HighestLag)
    {
        ind_lag_all <- lapply(1:N, function(x) {
            dayArgument_pers = if (!is.null(Day))
                Data[PersStart[x]:PersEnd[x], Day]
            else
                NULL
            out <-
                DetPredSubj(beep = Data[PersStart[x]:PersEnd[x], Beep],
                            day = dayArgument_pers,
                            MaxLag = lagRunner)
            return(out)
        })
        ind_lag_all <- do.call(c, ind_lag_all)
        stopifnot(length(which(is.na(
            ind_lag_all
        ))) == 0) # check that everything went well
        PredictableObs[[lagRunner]] = as.data.frame(cbind(ind_lag_all, Pers))
        stopifnot(all(PredictableObs[[lagRunner]][, 2] == as.character(Data[, ID]))) # a check if person indicators in PredictableObs are ordered correctly
    }


    ### Runners that depend on the number of predictable obs (and thus on the number of lags) ### ----
    Tni_NPred = vector("list", HighestLag)  # from 1:HighestLag but only LowestLag:HighestLag elements are filled
    for (lagRunner in LowestLag:HighestLag)
    {
        Tni_NPred[[lagRunner]] = rep(0, N)
        nPredObsPerPerson = PredictableObs[[lagRunner]][PredictableObs[[lagRunner]][, 1] == 1, 2]
        for (i in 1:N)
        {
            # determine the total number of predictable Observations per person (Tni_NPred) differs across lag numbers
            Tni_NPred[[lagRunner]][i] = length(which(nPredObsPerPerson == pers[i]))
            # Tni_NPred[[lagRunner]] # length of U per person,
            # With presample of first Lags-obs removed for every individual and all other obs that cannot be predicted removed,
        }

    }


    FunctionOutput = data.frame(matrix(NA, nrow = HighestLag, ncol = 1),
                                row.names = apply(as.matrix(1:HighestLag), 1, function(x) paste(c(x, "Lag"), collapse = " ")))
    colnames(FunctionOutput) = c("Total predictable observations")

    for (lagRunner in LowestLag:HighestLag)
    {
        FunctionOutput[lagRunner, ] = sum(Tni_NPred[[lagRunner]])
    }
    FunctionOutput = FunctionOutput[LowestLag:HighestLag, , drop = FALSE]
    #class(FunctionOutput) <- c("ClusterVARNumberPredictable", class(FunctionOutput))

    return(FunctionOutput)

}
