# jonashaslbeck@protonmail.com; Nov 13th, 2023

# --------------------------------------------------
# ---------- What is happening here? ---------------
# --------------------------------------------------

# Writing up model selection loop with AIC / BIC
# This should later be in a function

# Note: For now, this is just running over a sequence
# of number of clusters, but this should later go over
# a grid of cluster numbers and lag sizes


# --------------------------------------------------
# ---------- Loading Packages ----------------------
# --------------------------------------------------

library(RColorBrewer) # for coloors


# --------------------------------------------------
# ---------- Loading Data --------------------------
# --------------------------------------------------

# We just use the example data in the package, which
# I assume are synthetic K=2 data

head(ExampleData)
plot(ExampleData[1:200, 6]) # OK, so there are no gaps ...
N <- nrow(ExampleData)


# --------------------------------------------------
# ---------- Loop over K-sequence ------------------
# --------------------------------------------------

clusterSeq <- 1:4 # We search clusters 1, 2, 3, 4
m_IC <- matrix(NA, 4, 6)
colnames(m_IC) <- c("LL", "nPar", "BIC", "AIC", "IC", "SC")
# l_models <- list()

for(i in 1:4) {

  l_models[[i]] <- LCVARclust(Data = ExampleData,
                              yVars = 1:4,
                              Time = 6,
                              ID = 5,
                              Covariates = "equal-within-clusters",
                              Clusters = clusterSeq[i],
                              LowestLag = 1,
                              HighestLag = 1,
                              smallestClN = 3,
                              ICType = "HQ",
                              seme = 3,
                              Rand = 2, # number of random re-starts
                              Rational = TRUE,
                              SigmaIncrease = 10,
                              it = 25,
                              Conv = 1e-06,
                              xContinuous = 7,
                              xFactor = 8)

  m_IC[i,1] <- l_models[[i]]$BestSolutionsPerCluster[[1]]$last.loglik
  m_IC[i,2] <- l_models[[i]]$BestSolutionsPerCluster[[1]]$nPara

  m_IC[i,5] <- l_models[[i]]$BestSolutionsPerCluster[[1]]$IC
  m_IC[i,6] <- l_models[[i]]$BestSolutionsPerCluster[[1]]$SC



  # ModelSummaryLCVAR(l_models[[3]], Output = "HQ")

  print(i)

} # end for


# NOTE: This loop should be packaged in some function
#       And on top of this we'll have S3methods (like
#       print, summary, plot) that create summaries of
#       the model selection such as the below figure



# --------------------------------------------------
# ---------- Evaluate ------------------------------
# --------------------------------------------------

# nice colors
cols <- brewer.pal(4, "Set1")

# Compute BIC
m_IC[, 3] <- -2*m_IC[, 1] + m_IC[, 2] * log(N)
# Compute AIC
m_IC[, 4] <- -2*m_IC[, 1] + 2*m_IC[, 2]

# Normalize them all between 01
m_IC_norm_p <- apply(m_IC[, 3:6], 2, function(x) {
  y <- x - min(x)
  y / max(y)
}  )


# Plotting
par(mar=c(4.5,5,2,1.5))
plot.new()
# y_range <- range(m_IC[,3:4])
# alpha <- 0.05
# y_range[1] <- y_range[1] * (1-alpha)
# y_range[2] <- y_range[2] * (1+alpha)
plot.window(xlim=range(clusterSeq), ylim=c(0,.05))
axis(1, labels = clusterSeq, at=clusterSeq)
# axis(2, las=2, round(seq(y_range[1], y_range[2], length=8)))
axis(2, las=2)
title(xlab = "Number of Clusters")
title(ylab = "Normalized IC Value", line=4)

# BIC
lines(clusterSeq, m_IC_norm_p[, 1], col=cols[1], lwd=2)
K_sel_BIC <- which.min(m_IC_norm_p[, 1])
points(K_sel_BIC, m_IC_norm_p[K_sel_BIC, 1], cex=2, pch=21, col=cols[1])
# AIC
lines(clusterSeq, m_IC_norm_p[, 2], col=cols[2], lwd=2)
K_sel_AIC <- which.min(m_IC_norm_p[, 2])
points(K_sel_AIC, m_IC_norm_p[K_sel_AIC, 2], cex=2, pch=21, col=cols[2])
# IC
lines(clusterSeq, m_IC_norm_p[, 3], col=cols[3], lwd=2)
K_sel_IC <- which.min(m_IC_norm_p[, 3])
points(K_sel_IC, m_IC_norm_p[K_sel_IC, 3], cex=2, pch=21, col=cols[3])
# SC
lines(clusterSeq, m_IC_norm_p[, 4], col=cols[4], lwd=2)
K_sel_SC <- which.min(m_IC_norm_p[, 4])
points(K_sel_SC, m_IC_norm_p[K_sel_SC, 4], cex=2, pch=21, col=cols[4])


# Legend
legend("topright", c("BIC", "AIC", "IC", "SC"), text.col=cols, cex=1.5, bty="n")


# Note: when we search also lags, we need a different visualization, maybe with heat plots



