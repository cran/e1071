rmvnorm <- function(n, mu=rep(0, nrow(sigma)),
                      sigma=diag(length(mu))){

  if(nrow(sigma) != ncol(sigma)){
    stop("sigma must be a square matrix")
  }

  if(length(mu) != nrow(sigma)){
    stop("mu and sigma have non-conforming size")
  }
  
  sigsvd <- svd(sigma)
  retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
  retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
  retval <- sweep(retval, 2, mu, "+")
  retval
}


dmvnorm <- function(x, mu, sigma){

  if(is.vector(x)){
    x <- matrix(x, ncol=length(x))
  }

  if(missing(mu)){
    mu <- rep(0, length=ncol(x))
  }
  
  if(missing(sigma)){
    sigma <- diag(ncol(x))
  }

  if(length(x) != ncol(sigma)){
    stop("x and sigma have non-conforming size")
  }
  
  if(nrow(sigma) != ncol(sigma)){
    stop("sigma must be a square matrix")
  }
  if(length(mu) != nrow(sigma)){
    stop("mu and sigma have non-conforming size")
  }

  retval <- exp(-mahalanobis(x, center=mu, cov=sigma)/2)
  det <- prod(eigen(sigma, sym=TRUE)$values)
  retval<- retval / (sqrt(det) * sqrt(2*pi)^ncol(x))

  retval
}
  
