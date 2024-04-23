
plot.ClusterVAR <- function(x,
                            ...) {


  # ----- Fill in defaults ------
  args <- list(...)

  if(!args$show %in% c("GNC", "GNL", "specific", "specificDiff")) stop("Inadmissible input for argument 'show'. see ?plot.ClusterVAR")
  if(is.null(args$show)) stop("Select which plot should be shown (see ?plot.ClusterVAR).") else show <- args$show
  if(is.null(args$Number_of_Clusters)) Number_of_Clusters <- NULL else Number_of_Clusters <- args$Number_of_Clusters
  if(is.null(args$Number_of_Lags)) Number_of_Lags <- min(x$Call$Lags) else Number_of_Lags <- args$Number_of_Lags
  if(is.null(args$Model)) Model <- NULL else Model <- args$Model
  if(is.null(args$labels)) labels <- NULL else labels <- args$labels
  if(is.null(args$cex.axis)) cex.axis <- 0.8 else cex.axis <- args$cex.axis
  if(is.null(args$cex.val)) cex.val <- 0.7 else cex.val <- args$cex.val

  # ----- Get relevant data from summary method ------
  # (if relevant)

  if(show %in% c("GNC", "GNL")) {

    out_sum <- summary(object = x,
                       show = show,
                       Number_of_Clusters = Number_of_Clusters,
                       Number_of_Lags = Number_of_Lags)
    out_table <- out_sum$FunctionOutput
  }


  # ----- Get relevant data from coef method ------
  # (if relevant): plots specific VAR model

  if(show %in% c("specific", "specificDiff")) {

    l_coef <- coef(x, Model = Model)

  }


  # ----- Global Plotting Settings ------
  cols <- c("#E41A1C", "#377EB8", "#4DAF4A")
  old_par <- par() # save old graphic settings

  # ----- Plotting: Best-per-number-of-clusters ------

  # if(show == "BPC") {
  #
  #   K <- nrow(out_table)
  #   yrange <- range(out_table$`log-likelihood`)
  #
  #   # Plotting canvas
  #   par(mar=c(4.4, 5.5, 2, 1.2))
  #   plot.new()
  #   plot.window(xlim=c(1,K), ylim=yrange)
  #   axis(1)
  #   axis(2, round(seq(yrange[1], yrange[2], length=8)), las=2)
  #   title(xlab="Number of Clusters")
  #   title(ylab="Log-likelihood", line=4.5)
  #   title("Different Number of Clusters (each with best Lag model)", font.main=1)
  #
  #   # Plotting Data
  #   points(1:K, out_table$`log-likelihood`, col=cols[3], pch=19)
  #   lines(1:K, out_table$`log-likelihood`, col=cols[3])
  #
  #   # Legend
  #   legend("bottomright", legend=c("Log-likelihood"),
  #          lty=1, col=cols[3], text.col=cols[3],
  #          bty="n", cex=1.2, pch=c(19))
  #
  # } # end if


  # ----- Plotting: "Given-a-number-of-clusters" ------

  if(show == "GNC") {

    N_L <- nrow(out_table) # number of models

    # Cut x-axis labels out of rownames of summary table
    # labels <- out_table$Lags
    # labels <- substr(rownames(out_table), 6, nchar(rownames(out_table)))
    labels <- out_table$Lags[N_L:1]
    yrange <- range(c(out_table$HQ, out_table$SC))

    # Plotting canvas
    par(mar=c(4.4, 5.5, 2, 1.2))
    plot.new()
    plot.window(xlim=c(1,N_L), ylim=yrange)
    axis(1, labels = labels, at=1:N_L)
    axis(2, labels=round(seq(yrange[1], yrange[2], length=8), 4),
         at=seq(yrange[1], yrange[2], length=8), las=2)
    title(xlab="Lag-Combinations")
    title(ylab="Information Criterion", line=4.5)
    title(paste0("Lags Combinations for ", Number_of_Clusters, " Clusters"), font.main=1)

    # Data
    points(1:N_L, out_table$HQ[N_L:1], col=cols[1], pch=19)
    points(1:N_L, out_table$SC[N_L:1], col=cols[2], lty=2, pch=17)
    lines(1:N_L, out_table$HQ[N_L:1], col=cols[1])
    lines(1:N_L, out_table$SC[N_L:1], col=cols[2], lty=2)

    # Legend
    legend("topleft", legend=c("HQ", "SC"),
           lty=1:2, col=cols[1:2], text.col=cols[1:2],
           bty="n", cex=1.2, pch=c(19, 17))



  } # end if


  # ----- Plotting: "Given-a-number-of-lags" ------

  if(show == "GNL") {

    K <- nrow(out_table) # number of models
    labels <- out_table$Lags
    yrange <- range(c(out_table$BIC, out_table$ICL))

    # Plotting canvas
    par(mar=c(4.4, 5.5, 2, 1.2))
    plot.new()
    plot.window(xlim=c(1,K), ylim=yrange)
    axis(1, labels = labels, at=1:K)
    axis(2, round(seq(yrange[1], yrange[2], length=8)), las=2)
    title(xlab="Models with different #Clusters")
    title(ylab="Information Criterion", line=4.5)
    title(paste0("Different Number of Clusters (Fixed lag = ", Number_of_Lags, ")"), font.main=1)

    # Plotting Data
    points(1:K, out_table$ICL, col=cols[1], pch=19)
    points(1:K, out_table$BIC, col=cols[2], lty=2, pch=17)
    lines(1:K, out_table$ICL, col=cols[1])
    lines(1:K, out_table$BIC, col=cols[2], lty=2)

    # Legend
    legend("topright", legend=c("ICL", "BIC"),
           lty=1, col=cols[1:2], text.col=cols[1:2],
           bty="n", cex=1.2, pch=c(19, 17))

  } # end if


  # ----- Plotting: VAR Parameter Matrices ------

  if(show == "specific") {

    dims_phi <- dim(l_coef$VAR_coefficients)
    p <- dims_phi[1]
    if(dims_phi[1] != dims_phi[2]) stop("Currently only implemented for lag-1 models")
    K <- dims_phi[3]
    # browser()

    l_phi <- list()
    for(k in 1:K) l_phi[[k]] <-  l_coef$VAR_coefficients[, , k]

    # Decide on layout
    if(K == 1) par(mfrow=c(1,1))
    if(K == 2) par(mfrow=c(1,2))
    if(K %in% 3:4) par(mfrow=c(2,2))
    if(K %in% 5:9) par(mfrow=c(3,3))
    if(K > 10) {
      ldim <- ceiling(sqrt(K))
      par(mfrow=c(ldim,ldim))
    }

    # browser()

    # Loop over clusters & plot
    for(k in 1:K) plotHeat(phi = l_phi[[k]],
                           k = k,
                           main = paste0("Cluster ", k),
                           labels = labels,
                           cex.axis = cex.axis,
                           cex.val = cex.val)

  } # end if

  # ----- Plotting: Cluster-Differences in VAR Parameter Matrices ------

  if(show == "specificDiff") {

    # Get parameters
    dims_phi <- dim(l_coef$VAR_coefficients)
    p <- dims_phi[1]
    if(dims_phi[1] != dims_phi[2]) stop("Currently only implemented for lag-1 models")
    K <- dims_phi[3]
    l_phi <- list()
    for(k in 1:K) l_phi[[k]] <-  l_coef$VAR_coefficients[, , k]


    # Setup Layout matrix
    lmat <- matrix((2*(K-1)+1):((2*(K-1)) + (K-1)^2), K-1, K-1, byrow=TRUE)
    lmat <- rbind(1:(K-1), lmat)
    lmat <- cbind(c(0, K:(2*(K-1))), lmat)
    layout(mat=lmat,
           widths = c(0.2, rep(1, K-1)),
           heights =  c(0.2, rep(1, K-1)))

    # Plot Labels
    plotLabel <- function(x, srt=0, col="black",
                          xpos=.6, ypos=.6, cex=1.4) {
      par(mar=rep(0, 4))
      plot.new()
      plot.window(xlim=c(0,1), ylim=c(0,1))
      text(xpos, ypos, x, srt=srt, cex=cex, col=col)
    }
    for(k in 2:K) plotLabel(paste0("Cluster ", k), cex=1.5)
    for(k in 1:(K-1)) plotLabel(paste0("Cluster ", k), cex=1.5, srt=90)

    # Plot Data
    for(k1 in 1:(K-1)) {
      for(k2 in 2:K) {

        if(k1==k2) {
          plot.new()
          plot.window(xlim=c(0,1), ylim=c(0,1))
        } else {
          phi_diff <- l_phi[[k1]] - l_phi[[k2]]
          par(mar=c(2.5,2.5,2,1))
          plotHeat(phi = phi_diff, k = k, main = paste0("Difference: Cluster ", k1, " - Cluster ", k2))
        }


      }
    } # end for Ks



  } # end if


  # Set graphic settings back to initial
  par(mfrow=old_par$mfrow, mar=old_par$mar)


} # eoF















