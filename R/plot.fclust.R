plot.fclust <- function(clobj, x, centers=TRUE, initcenters=FALSE,
                         color=rainbow(clobj$learning$ncenters),...){
  
  x <- as.matrix(x)

  
  if(dim(x)[2]>2){
    pairs(x, col=color[cl$cluster], ...)
  }
  else{
    plot(x, col=color[cl$cluster], ...)
    if(centers)
      points(cl$centers, pch=4,col=color,cex=2)
    if(initcenters)
      points(clobj$learning$initcenters, pch=3,col=color,cex=2)
  }
}

