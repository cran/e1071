kurtosis <- function (x, na.rm = FALSE)
{
  if (na.rm) 
    x <- x[!is.na(x)]
  sum((x-mean(x))^4)/(length(x)*var(x)^2) - 3
}


