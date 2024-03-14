

plotHeat <- function(phi, k, main) {

  p <- ncol(phi)

  # -- Make color gradient --
  color.gradient <- function(x, colors=c("#E41A1C","#377EB8"), colsteps=201) {
    return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
  }
  x <- 1:201
  grad <- color.gradient(x)
  # grad <- scales::alpha(grad, .75)

  # -- Plotting --
  # Make canvas
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









