print.fclust <- function (x, ...)
  {
    clobj <- x
    if (!is.null(clobj$iter))
      cat("\n                            Clustering on Training Set\n\n\n")
    else
      cat("\n                              Clustering on Test Set\n\n\n")
    
    cat("Number  of  Clusters: ", clobj$ncenters, "\n")
    cat("Sizes   of  Clusters: ", clobj$size, "\n\n")
    cat("Centers of  Clusters: ", clobj$centers, "\n\n")
    
    
#    cat("Learning Parameters:",clobj$call,"\n\n")
    
    if (clobj$iter < clobj$call$iter.max)
      cat("Algorithm converged after", clobj$iter, "iterations.\n")
    else
      cat("Algorithm did not converge after", clobj$iter, "iterations.\n")
  }


