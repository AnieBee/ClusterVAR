

plotHeat <- function(phi,
                     k,
                     main,
                     labels=NULL,
                     las.x=1,
                     cex.axis=0.8,
                     cex.val=0.7,
                     mar = c(2.5,2.5,2,1)) {


  # -- Aux Variables --
  p <- ncol(phi)

  # Fill in default labels
  if(is.null(labels)) labels <-  paste0("Y", 1:p)


  # -- Recover graphics settings --
  # Jonas Oct 9th: DON'T DO THIS HERE, because otherwise the layouts get reset when we repeatedly call this
  # oldpar <- par(no.readonly = TRUE) # code line i
  # on.exit(par(oldpar)) # code line i + 1

  # -- Make color gradient --
  color.gradient <- function(x, colors=c("#E41A1C", "white", "#377EB8"), colsteps=201) {
    return( grDevices::colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
  }
  x <- 1:201
  grad <- color.gradient(x)

  # Make canvas
  par(mar=mar)
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0, 1))

  # Auxiliary plotting variables
  sfm <- 1/(p*2)
  seq_mp_x <- seq(0, 1, length=p+1)[-(p+1)] + sfm

  # Plot Axes & Axis labels
  labels_tm1 <- sapply(labels, function(label) bquote(.(label)[t-1]), simplify = FALSE)
  labels_t1 <- sapply(labels, function(label) bquote(.(label)[t]), simplify = FALSE)
  axis(1, labels = do.call(expression, labels_tm1), at=seq_mp_x, cex.axis=cex.axis)
  axis(2, labels = do.call(expression, labels_t1[p:1]), at=seq_mp_x, las=2, cex.axis=cex.axis)
  title(main, font.main=1)

  # Plot Data
  for(i in 1:p) {
    for(j in 1:p) {

      # Get color
      phi_ij <- phi[p:1, ][j, i]
      if(phi_ij < -1) {
        col_ij <- grad[1]
      } else if(phi_ij > 1 ) {
        col_ij <- grad[201]
      } else {
        col_ij <- grad[phi[p:1, ][j, i] * 100 + 101]
      }

      # Plot box
      rect(xleft = seq_mp_x[i]-sfm,
           ybottom = seq_mp_x[j]-sfm,
           xright = seq_mp_x[i]+sfm,
           ytop = seq_mp_x[j]+sfm,
           col = col_ij)

      # Plot text
      text(seq_mp_x[i], seq_mp_x[j], round(phi_ij , 2), cex=cex.val, col="black")
    }
  }


} # eoF





