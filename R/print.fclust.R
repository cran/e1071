print.fclust <- function (x, ...)
  {
    clobj <- x
    if (!is.null(clobj$iter))
      cat("\n                            Clustering on Training Set\n\n\n")
    else
      cat("\n                              Clustering on Test Set\n\n\n")
    
    cat("Number  of  Clusters: ", dim(clobj$centers)[1], "\n")
    cat("Sizes   of  Clusters: ", clobj$size, "\n\n")
#    cat("Centers of  Clusters: ", "\n")
    #matrix(clobj$centers, ncol=dim(clobj$ce)[2]), "\n\n")
#    clobj$centers
    
#    cat("Learning Parameters:",clobj$call,"\n\n")
    
#    if (clobj$iter < clobj$call$iter.max)
#      cat("Algorithm converged after", clobj$iter, "iterations.\n")
#    else
#      cat("Algorithm did not converge after", clobj$iter, "iterations.\n")
  }


