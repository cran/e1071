read.matrix.csr <- function(file, fac = TRUE, ncol = NULL)
{
  if (!require(SparseM)) 
    stop("Need `SparseM' package!")
  l <- strsplit(readLines(file(file)), "[ ]+")

  ## y-values
  y <- sapply(l, function(x) x[1])

  ## x-values
  rja <- do.call("rbind", lapply(l, function(x) do.call("rbind", strsplit(x[-1], ":"))))
  ja <- as.integer(rja[,1])
  ia <- cumsum(c(1, sapply(l, length) - 1))

  max.ja <- max(ja)
  dimension <- c(length(l), if (is.null(ncol)) max.ja else max(ncol, max.ja))
  x = new("matrix.csr", as.numeric(rja[,2]), ja, as.integer(ia), as.integer(dimension))
  if (length(y)) 
    list(x = x, y = if (fac) as.factor(y) else as.numeric(y))
  else x
}

## old version: slow, but uses less memory
# read.matrix.csr <- function (file, fac = TRUE, ncol = NULL) 
# {
#     library(methods)
#     if (!require(SparseM)) 
#         stop("Need `SparseM' package!")
#     con <- file(file)
#     open(con)
#     y <- vector()
#     ia <- 1
#     ra <- ja <- c()
#     i <- 1
#     maxcol <- 1
#     while (isOpen(con) & length(buf <- readLines(con, 1)) > 0) {
#         s <- strsplit(buf, "[ ]+", extended = TRUE)[[1]]

#         ## y
#         if (length(grep(":", s[1])) == 0) {
#             y[i] <- if (fac) 
#                 s[1]
#             else as.numeric(s[1])
#             s <- s[-1]
#         }

#         ## x
#         if (length(s)) {
#             tmp <- do.call("rbind", strsplit(s, ":"))
#             ra <- c(ra, as.numeric(tmp[, 2]))
#             ja <- c(ja, as.numeric(tmp[, 1]))
#         }
        
#         i <- i + 1
#         ia[i] <- ia[i - 1] + length(s)
#     }
#     dimension <- c(i - 1, if (is.null(ncol)) max(ja) else max(ncol, ja))
#     x = new("matrix.csr", ra, as.integer(ja), as.integer(ia), as.integer(dimension))
#     if (length(y)) 
#         list(x = x, y = if (fac) as.factor(y) else y)
#     else x
# }

write.matrix.csr <- function (x, file="out.dat", y=NULL) {
  library(methods)
  if (!is.null(y) & (length(y) != nrow(x)))
    stop(paste("Length of y (=", length(y),
                 ") does not match number of rows of x (=",
                 nrow(x), ")!", sep=""))
  sink(file)
  for (i in 1:nrow(x)) {
    if (!is.null(y)) cat (y[i],"")
    for (j in x@ia[i]:(x@ia[i+1] - 1))
      cat(x@ja[j], ":", x@ra[j], " ", sep="")
    cat("\n")
  }
  sink()
}

na.fail.matrix.csr <- function(object, ...) {
  library(methods)
  if (any(is.na(object@ra)))
    stop("missing values in object") else return(object)
}






