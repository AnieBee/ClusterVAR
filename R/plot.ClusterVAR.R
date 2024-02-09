


plot.ClusterVAR <- function(object,
                            show = "Best-per-number-of-clusters",
                            Number_of_Clusters = NULL) {

  # ----- Get relevant data from summary method ------

  out_sum <- summary(LCVARresult = object,
                     show = show,
                     Number_of_Clusters = Number_of_Clusters)
  out_table <- out_sum$FunctionOutput


  # ----- Global Plotting Settings ------
  cols <- c("#E41A1C", "#377EB8", "#4DAF4A")

  # ----- Plotting: Best-per-number-of-clusters ------

  if(show == "Best-per-number-of-clusters") {

    K <- nrow(out_table)
    yrange <- range(c(out_table$ICL, out_table$BIC))

    # Plotting canvas
    par(mar=c(4.4, 5.5, 2, 1.2))
    plot.new()
    plot.window(xlim=c(1,K), ylim=yrange)
    axis(1)
    axis(2, round(seq(yrange[1], yrange[2], length=8)), las=2)
    title(xlab="Number of Clusters")
    title(ylab="Information Criterion", line=4.5)
    title("Models with Best Lag for each Cluster", font.main=1)

    # Plotting Data
    lines(1:K, out_table$ICL, col=cols[1])
    lines(1:K, out_table$BIC, col=cols[1], lty=2)

    # Legend
    legend("topright", legend=c("ICL", "BIC"),
           lty=1:2, col=cols[1:2], text.col=cols[1:2],
           bty="n", cex=1.2)

  } # end if


  # ----- Plotting: "Given-a-number-of-clusters" ------

  if(show == "Given-a-number-of-clusters") {

    N_L <- nrow(out_table) # number of models

    # TODO: Order table by complexity
    labels <- substr(rownames(out_table), 6, nchar(rownames(out_table)))
    yrange <- range(c(out_table$HQ, out_table$SC))

    # Get labels ..

    # Plotting canvas
    par(mar=c(4.4, 5.5, 2, 1.2))
    plot.new()
    plot.window(xlim=c(1,N_L), ylim=yrange)
    axis(1, labels = labels, at=1:N_L)
    axis(2, labels=round(seq(yrange[1], yrange[2], length=8), 4),
         at=seq(yrange[1], yrange[2], length=8), las=2)
    title(xlab="Lags Combinations")
    title(ylab="Information Criterion", line=4.5)
    title(paste0("Lags Combinations for ", Number_of_Clusters, " Clusters"), font.main=1)

    # Data
    lines(1:N_L, out_table$HQ, col=cols[1])
    lines(1:N_L, out_table$SC, col=cols[2], lty=2)

    # Legend
    legend("topleft", legend=c("HQ", "SC"),
           lty=1:2, col=cols[1:2], text.col=cols[1:2],
           bty="n", cex=1.2)



  } # end if


  # ----- Plotting: XXX ------



} # eoF
