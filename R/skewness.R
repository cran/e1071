skewness <- function (x, na.rm = FALSE)
{
  if (na.rm) 
    x <- x[!is.na(x)]
  sum((x-mean(x))^3)/(length(x)*sd(x)^3)
}

