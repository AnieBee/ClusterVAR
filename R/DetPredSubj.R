# jonashaslbeck@protonmail.com; Feb 16th, 2024

# --------------------------------------------------------
# ---------- What is happening here? ---------------------
# --------------------------------------------------------

# Make function to determine whether a given time point can
# be predicted given that we git a VAR(s) model

# --------------------------------------------------------
# ---------- Function to determine valid time points -----
# --------------------------------------------------------

# Input:
# 1) n x p data matrix, including day number, beep number, subj id, modeled variables
# 2) lag order s
# Output
# 1) n x 1 vector, indicating whether each time point can be predicted

DetPredSubj <- function(beep, day, MaxLag=1) {
    
    Nt <- length(beep)
    ind_pred <- rep(NA, Nt)
    ind_pred[1:MaxLag] <- FALSE # First MaxLag can never be predicted
    for(i in (MaxLag+1):Nt) { 
        # Check lags
        v_check1 <- rep(NA, MaxLag)
        for(j in 1:MaxLag) v_check1[j] <- (beep[i] - beep[i-j]) == j
        check1 <- all(v_check1)
        
        if(!is.null(day)){
            # Same day?
            check2 <- day[i] == day[i-1]
            # Make time vector
            ind_pred[i] <- ifelse(check1 & check2, 1, 0)
        }else{
            ind_pred[i] <- ifelse(check1, 1, 0)
        }
    } # end for
    return(ind_pred)
} # eoF


