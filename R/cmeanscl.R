cmeanscl <- function (x, centers, iter.max = 100, verbose = FALSE,
                    dist = "euclidean", method = "cmeans",
                    m=2, rate.par = NULL) 
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
  
  dist <- pmatch(dist, c("euclidean", "manhattan"))
  if (is.na(dist)) 
    stop("invalid distance")
  if (dist == -1) 
    stop("ambiguous distance")
  
  method <- pmatch(method, c("cmeans", "ufcl"))
  if (is.na(method)) 
    stop("invalid clustering method")
  if (method == -1) 
    stop("ambiguous clustering method")
  
  if (method == 2) {
    if (missing(rate.par)) {
      rate.par <- 0.3
    }
  }
            
  initcenters <- centers
  ## dist <- matrix(0, xrows, ncenters)
  ## necessary for empty clusters
  pos <- as.factor(1:ncenters)
  rownames(centers) <- pos
  iter <- integer(1)
  
  if (method == 1) {
    retval <- .C("cmeans", xrows = as.integer(xrows),
                 xcols = as.integer(xcols), 
                 x = as.double(x), ncenters = as.integer(ncenters), 
                 centers = as.double(centers), 
                 iter.max = as.integer(iter.max), iter = as.integer(iter), 
                 verbose = as.integer(verbose), dist = as.integer(dist-1), 
                 U=double(xrows*ncenters), UANT=double(xrows*ncenters),
                 m=as.double(m), ermin=double(1))
  }
  else if (method == 2) {
    retval <- .C("ufcl", xrows = as.integer(xrows),
                 xcols = as.integer(xcols), 
                 x = as.double(x), ncenters = as.integer(ncenters), 
                 centers = as.double(centers), 
                 iter.max = as.integer(iter.max), iter = as.integer(iter), 
                 verbose = as.integer(verbose), dist = as.integer(dist-1), 
                 U=double(xrows*ncenters), m=as.double(m),
                 rate.par = as.double(rate.par), ermin=double(1))
  }
  
  centers <- matrix(retval$centers, ncol = xcols, dimnames = dimnames(initcenters))
  
  U <- retval$U
  U <- matrix(U, ncol=ncenters)
#  clusterU <- max.col(U)
   clusterU <- apply(U,1,which.max)
  clusterU <- clusterU[order(perm)]
  U <- U[order(perm),]
  
  clustersize <- as.integer(table(clusterU))
  
  retval <- list(centers = centers, cluster = clusterU,
                 size = clustersize, dist = dist,
                 m=retval$m, member=U, withinsd = retval$ermin,
                 learning = list(ncenters = ncenters,
                   initcenters = initcenters, iter = retval$iter - 1,
                   rate.par = rate.par), call = match.call())
  
  class(retval) <- c("fclust")
  return(retval)
}





