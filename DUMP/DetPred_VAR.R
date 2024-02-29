# jonashaslbeck@protonmail.com; Feb 16th, 2024

# --------------------------------------------------------
# ---------- What is happening here? ---------------------
# --------------------------------------------------------

# Make function to determine whether a given time point can
# be predicted given that we git a VAR(s) model

# Note: some datasets have already missing time points removed
# while others have all the rows in there, but then responses are
# missing; in the code below I remove missing values first, and
# thereby we end up with the same datastructre no matter what
# format we look at
# Of course, in the package we want to have input checks around this
# so we are sure the data looks like what we want it to look like


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
    
    # Same day?
    check2 <- day[i] == day[i-1]
    
    # Make time vector
    ind_pred[i] <- ifelse(check1 & check2, TRUE, FALSE)
    
  } # end for
  
  return(ind_pred)
  
} # eoF



# --------------------------------------------------------
# ---------- Showcase ------------------------------------
# --------------------------------------------------------

# ------ Load Data of Fried et al. (2021) ------

data <- readRDS("data_Fried2021.RDS")


# ------ Showcase single subject ------

# Take the first subject
u_subj <- unique(data$subj_id)
data_i <- data[data$subj_id == u_subj[1], ]
nrow(data_i)

# Get rid of NAs
data_i_noNA <- na.omit(data_i) # Of course this requires that there are no NAs in some nuisance columns!
nrow(data_i_noNA)

# Apply function: Lag=1
ind_lag1 <- DetPredSubj(beep=data_i_noNA$beepvar, 
                        day=data_i_noNA$dayvar, 
                        MaxLag=1)

length(ind_lag1)

# Let's have a look
inspect <- cbind(data_i_noNA$dayvar, data_i_noNA$beepvar, ind_lag1)
colnames(inspect) <- c("day", "beep", "ind")
inspect
# Some examples: 
# 1) row 1 can never be predicted
# 2) row 5 cannot be predicted because of the day change
# 3) row 9 cannot be predicted, because, while we are on the same day (3), we only have measurements 1 and 4; 
#.   to predict 4 with a lag 1 model we would need measurement 3 though, which is missing


# Apply function: Lag=1
ind_lag2 <- DetPredSubj(beep=data_i_noNA$beepvar, 
                        day=data_i_noNA$dayvar, 
                        MaxLag=2)

length(ind_lag2)
inspect <- cbind(data_i_noNA$dayvar, data_i_noNA$beepvar, ind_lag2)
colnames(inspect) <- c("day", "beep", "ind")
inspect
# Some examples:
# 1) row 1 and 2 can never be predicted in Lag 2 model
# 2) row 17 could be predicted in Lag(1) but not Lag(2)

# ------ Apply to whole dataset ------

# Omit missing
data_noNA <- na.omit(data)
# Unique participants
u_subj <- unique(data$subj_id)
# Loop with sapply
ind_lag_1_all <- sapply(1:length(u_subj), function(x) {
  data_i <- data_noNA[data_noNA$subj_id==u_subj[x], ]
  out <- DetPredSubj(beep = data_i$beepvar, 
                     day = data_i$dayvar, 
                     MaxLag = 1)
  return(out)
})
v_ind_lag_1_all <- do.call(c, ind_lag_1_all)

nrow(data_noNA)
length(v_ind_lag_1_all) # match!


















