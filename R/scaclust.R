scaclust <- function (x, centers, iter.max = 100, verbose = FALSE,
                    method = "ad", theta= NULL) 
{
  xrows <- dim(x)[1]
  xcols <- dim(x)[2]
  xold <- x
  perm <- sample(xrows)
  x <- x[perm, ]
  ## initial values are given
  if (is.matrix(centers)) 
    ncenters <- dim(centers)[1]
  else {
   ## take centers random vectors as initial values
    ncenters <- centers
    centers <- x[rank(runif(xrows))[1:ncenters], ]+0.001
  }
##  method <- pmatch(method, c("ad", "mtv", "sand","ml"))
  method <- pmatch(method, c("ad", "mtv", "sand","mlm"))
  if (is.na(method)) 
    stop("invalid clustering method")
  if (method == -1) 
    stop("ambiguous clustering method")
  
  if (method == 1) {
    beta <- 1/xcols
    taf <- 0 }
  if (method == 2) {
    beta <- 0.5
    taf <- xcols/2 }
  if (method == 3) {
    beta <- 1/xcols
    taf <- 1 }
  if (method == 4){
    beta <- 0.0
    taf <- -1 }

  
  ##initialize theta
 ## if (method != 4){
  if (missing(theta))
    theta <- rep(1.0,ncenters)
  else
    theta <- as.double(theta)
  
  ##}
  
  
  initcenters <- centers
  ## dist <- matrix(0, xrows, ncenters)
  ## necessary for empty clusters
  pos <- as.factor(1:ncenters)
  rownames(centers) <- pos
  iter <- integer(1)
  
##  if ((method == 1) || (method == 2) || (method == 3)){
    retval <- .C("common", xrows = as.integer(xrows),
                 xcols = as.integer(xcols), 
                 x = as.double(x), ncenters = as.integer(ncenters), 
                 centers = as.double(centers), 
                 itermax = as.integer(iter.max), iter = as.integer(iter), 
                 verbose = as.integer(verbose), U=double(xrows*ncenters),
                 UANT=double(xrows*ncenters), beta=as.double(beta),
                 taf=as.double(taf), theta=as.double(theta),ermin=double(1))

  centers <- matrix(retval$centers, ncol = xcols, dimnames = dimnames(initcenters))
  
  U <- retval$U
  U <- matrix(U, ncol=ncenters)
  ##  clusterU <- max.col(U)
  clusterU <- apply(U,1,which.max)
  clusterU <- clusterU[order(perm)]
  U <- U[order(perm),]
  
  clustersize <- as.integer(table(clusterU))
  
  retval <- list(centers = centers, cluster = clusterU,
                 size = clustersize, member=U, error = retval$ermin,
                 learning = list(ncenters = ncenters,
                   initcenters = initcenters, iter = retval$iter - 1,
                   theta = theta), call = match.call())
  
  class(retval) <- c("fclust")
  return(retval)
}



