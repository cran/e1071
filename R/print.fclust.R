print.fclust <- function (clobj)
  {
    
    if (!is.null(clobj$learning$iter))
      cat("\n                            Clustering on Training Set\n\n\n")
    else
      cat("\n                              Clustering on Test Set\n\n\n")
    
    cat("Number of Clusters: ", clobj$learning$ncenters, "\n")
    cat("Sizes  of Clusters: ", clobj$size, "\n\n")
    
  }

