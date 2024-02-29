# jonashaslbeck@gmail.com

# ---------------------------------------------------
# ---------- What is happening here? ----------------
# ---------------------------------------------------

# Minimal example for how to use foreach


# ---------------------------------------------------
# ---------- Load packages --------------------------
# ---------------------------------------------------

# Of course, in the package this has to be loaded via Namespace
library(foreach)
library(parallel)
library(doParallel)


# ---------------------------------------------------
# ---------- Example --------------------------------
# ---------------------------------------------------

nIter <- 6

# ----- for-loop ------
timer <- proc.time()[3]
v_out <- rep(NA, nIter)
for(i in 1:nIter) {
  set.seed(i)
  data <- rnorm(5*10^7)
  v_out[i] <- mean(data)
}
proc.time()[3] - timer


# ----- parallelized foreach-loop ------
timer <- proc.time()[3]
cluster <- 6
cl <- makeCluster(cluster, outfile="")
registerDoParallel(cl)
# Below: For this examples I don't need to send objects or packages to the nodes, so I left this
# commented out; but this shows how to do it

v_out2 <- foreach(i = 1:nIter
                  # .packages = c("bgms", "IsingFit"),
                  # .export = c("n_seq", "v_n_edges_present", "n_methods", "n_dvars", "GenData_RndG", "EstimationFunction"),
                  ) %dopar% {

  set.seed(i)
  data <- rnorm(5*10^7)
  return(mean(data))

} # eoF
stopCluster(cl)
proc.time()[3] - timer

# Note: here the computational cost of the random number generator is relatively small compared
# to the overhead of parallelization; if the latter becomes larger, the parallelization advantage
# should grow towards 6:1 (in this case of 6 cores)

# Same results:
v_out
unlist(v_out2)










