my.biplot.pcoa<-function (x, Y = NULL, plot.axes = c(1, 2), dir.axis1 = 1, dir.axis2 = 1, rn = NULL,mar= c(4, 4, 1, 0.2), sameAxis=TRUE, ...) {
  k <- ncol(x$vectors)
  if (k < 2) 
    stop("There is a single eigenvalue. No plot can be produced.")
  if (k < plot.axes[1]) 
    stop("Axis", plot.axes[1], "does not exist.")
  if (k < plot.axes[2]) 
    stop("Axis", plot.axes[2], "does not exist.")
  if (!is.null(rn)) 
    rownames(x$vectors) <- rn
  labels = colnames(x$vectors[, plot.axes])
  diag.dir <- diag(c(dir.axis1, dir.axis2))
  x$vectors[, plot.axes] <- x$vectors[, plot.axes] %*% diag.dir
  if (is.null(Y)) {
    limits <- apply(x$vectors[, plot.axes], 2, range)
    ran.x <- limits[2, 1] - limits[1, 1]
    ran.y <- limits[2, 2] - limits[1, 2]
    xlim <- c((limits[1, 1] - ran.x/10), (limits[2, 1] + 
        ran.x/5))
    ylim <- c((limits[1, 2] - ran.y/10), (limits[2, 2] + 
        ran.y/10))
    par(mai = c(1, 1, 1, 0.5))
    plot(x$vectors[, plot.axes], xlab = labels[1], ylab = labels[2], 
      xlim = xlim, ylim = ylim, asp = 1)
    text(x$vectors[, plot.axes], labels = rownames(x$vectors), 
      pos = 4, cex = 1, offset = 0.5)
    title(main = "PCoA ordination", line = 2.5)
  }
  else {
    n <- nrow(Y)
    points.stand <- scale(x$vectors[, plot.axes])
    S <- cov(Y, points.stand)
    U <- S %*% diag((x$values$Eigenvalues[plot.axes]/(n - 
          1))^(-0.5))
    colnames(U) <- colnames(x$vectors[, plot.axes])
    par(mar = mar)
    myBiplot(x$vectors[, plot.axes], U, xlab = labels[1], ylab = labels[2],sameAxis=sameAxis,...)
  }
  invisible(x$vectors[, plot.axes])
}

myBiplot.pca<-myBiplot<-function(pcaPoints,pcaArrows,choice=c(1,2),cor=FALSE,prcomp=FALSE,arrowsFilter=NULL,limScale=1,sameAxis=TRUE,...){
  if(sameAxis){
    xlim<-range(pcaPoints) #find limits for x and y directions
    ylim<-xlim
  }else{
    xlim<-range(pcaPoints[,1])
    ylim<-range(pcaPoints[,2])
  }
	arrowScale<-max(range(pcaArrows[,1])/xlim,range(pcaArrows[,2])/ylim)*limScale #find a scaling for the arrows limits
	if(!is.null(arrowsFilter)){pcaArrows<-pcaArrows[sqrt(pcaArrows[,1]^2+pcaArrows[,2]^2)>arrowsFilter,,drop=FALSE]}
  if(nrow(pcaArrows)>0){
    plot(pcaArrows,type="n",xaxt='n',yaxt='n',xlab='',ylab='',xlim=arrowScale*xlim,ylim=arrowScale*ylim) #set the axes for easy arrow drawing (messes up any additional plotting on score axes)
    arrows(0,0,pcaArrows[,1]*.8,pcaArrows[,2]*.8,length=.1) #draw arrows
    arrowText<-dimnames(pcaArrows)[[1]] #get rownames of loadings for labels
    if(is.null(arrowText)) arrowText<-paste('Var',1:nrow(pcaArrows)) #if no rownames just label with "Var X"
    text(pcaArrows,arrowText,cex=.9) #label the arrows
    par(new=TRUE) #plot the next plot directly on top of the current one
  }
	plot(pcaPoints,xlim=xlim,ylim=ylim,...)
}

