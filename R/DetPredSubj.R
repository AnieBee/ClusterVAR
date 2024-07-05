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

DetPredSubj <- function(beep, day, MaxLag = 1) {
  Nt <- length(beep)
  ind_pred <- rep(NA, Nt)
  ind_pred[1:MaxLag] <- FALSE # First MaxLag can never be predicted

  for(i in (MaxLag+1):Nt) {
    check1 <- check2 <- check3 <- FALSE
    if(!is.null(day)){
      check1 <- (day[i-1] == day[i]) & ((beep[i] - beep[i-1]) == 1)
      if(MaxLag > 1) check2 <- (day[i-2] == day[i]) & ((beep[i] - beep[i-2]) == 2)
      if(MaxLag > 2) check3 <- (day[i-3] == day[i]) & ((beep[i] - beep[i-3]) == 3)
    } else{
      check1 <- ((beep[i] - beep[i-1]) == 1)
      if(MaxLag > 1) check2 <-((beep[i] - beep[i-2]) == 2)
      if(MaxLag > 2) check3 <- ((beep[i] - beep[i-3]) == 3)
    }

    if(MaxLag==1) ind_pred[i] = as.numeric(check1)
    if(MaxLag==2) ind_pred[i] = as.numeric(check1 & check2)
    if(MaxLag==3) ind_pred[i] = as.numeric(check1 & check2 & check3)
  }

  return(ind_pred)
} # eoF


