

plotHeat <- function(phi, k, main) {

  p <- ncol(phi)

  # browser()

  # -- Make color gradient --
  color.gradient <- function(x, colors=c("#E41A1C","#377EB8"), colsteps=201) {
    return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
  }
  x <- 1:201
  grad <- color.gradient(x)
  # grad <- scales::alpha(grad, .75)

  # # -- Make Legend --
  # # Create a sequence from -1 to 1
  # f_colLeg <- function() {
  #   par(mar=c(2.5,1,2,2.5))
  #   values <- seq(-1, 1, length.out = max(x))
  #   # Generate the colors for the gradient
  #   colors <- color.gradient(values)
  #   # Create an empty plot
  #   plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 1), ylim = c(-1, 1), axes = FALSE)
  #   # Add the color gradient as a vertical legend
  #   for (i in 1:(length(values) - 1)) rect(0, values[i], 1, values[i+1], col = colors[i], border = NA)
  #   x_axis <- seq(-1, 1, length=5)
  #   axis(4, at=x_axis, labels=x_axis, las=2, cex.axis=0.8)
  #   # for (i in 1:length(x_axis)) mtext(text = x_axis[i], side = 4, at = x_axis[i], las = 2, cex = .75, adj = 0, line=3) # Increase adj to move text left
  # }
  #
  # # -- Plotting --
  # # Make Layout
  # lmat <- matrix(2:1, nrow=1)
  # layout(mat = lmat, widths = c(1,.2))
  #
  #
  # # Plot Legend
  # f_colLeg()

  # Make canvas
  par(mar=c(2.5,2.5,2,1))
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0, 1))

  # Auxiliary plotting variables
  sfm <- 1/(p*2)
  seq_mp_x <- seq(0, 1, length=p+1)[-(p+1)] + sfm

  # Plot Axes
  axis(1, labels = paste0("X", 1:p, "(t-1)"), at=seq_mp_x, cex.axis=0.8)
  axis(2, labels = paste0("X", p:1, "(t)"), at=seq_mp_x, las=2, cex.axis=0.8)
  title(main, font.main=1)

  # Plot Data
  for(i in 1:p) {
    for(j in 1:p) {
      rect(xleft = seq_mp_x[i]-sfm,
           ybottom = seq_mp_x[j]-sfm,
           xright = seq_mp_x[i]+sfm,
           ytop = seq_mp_x[j]+sfm,
           col = grad[phi[p:1, ][j, i] * 100 + 101])
      text(seq_mp_x[i], seq_mp_x[j], round(phi[p:1, ][j,i] , 2), cex=.7, col="black")
    }
  }


} # eoF









