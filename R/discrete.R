rdiscrete <- function (n, probs, values = 1:length(probs), method="inverse",
                       aliasmatrix = NULL)
{
    
    if (length(probs) != length(values))
        stop("rdiscrete: probs and values must have the same length.")
    if (sum(probs < 0) > 0)
        stop("rdiscrete: probs must not contain negative values.")
    
    if (n == 1)
        return (values[sum(runif(1) > p) + 1])
    else
    {
        method <- pmatch(method, c("inverse", "alias"))
        if (is.na(method))
            stop("rdiscrete: unknown method.")
        if (method == 1)
        {
            p <- cumsum(probs)/sum(probs)
            l <- length(probs)
            m <- numeric(n)
            a <- runif(n)
            for (i in 1:n)
                m[i] <- sum(a[i] > p)
            return(values[m + 1])
        }
        else
        {
            if (missing(aliasmatrix))
                aliasmatrix <- aliasmat(probs)
            
            x <- sample(1:nrow(aliasmatrix), n, replace=TRUE)
            y <- runif(n)
            
            retval <- rep(0, length=n)
    
            eins <- (y <= aliasmatrix[x,1])
    
            retval[eins] <- values[aliasmatrix[x[eins],2]]
            retval[!eins] <- values[aliasmatrix[x[!eins],3]]
            
            return(retval)
        }
    }
}

aliasmat <- function(p)
{
    p <- p / sum(p)
    q <- p * (pn <- length(p))
    r <- matrix(0, nrow=pn, ncol=3)

    eps <- .Machine$double.eps
    while (sum(!is.na(q))>1)
    {
        qklein <- min((1:pn)[q<=1+eps],na.rm=TRUE)
        qgross <- max((1:pn)[q>=1-eps],na.rm=TRUE)
        r[qklein,] <- c(q[qklein],qklein,qgross)
        q[qgross] <- q[qgross] + q[qklein] - 1
        q[qklein] <- NA
    }
    qmittel <- (1:pn)[!is.na(q)]
    r[qmittel,] <-  c(q[qmittel],qmittel,qmittel)

    return(r)
}

    
aliasmat2prob <- function(aliasmatrix)
{
    p <- rep(0, length = length(unique(aliasmatrix[,2:3])))
    names(p) <- (pnames <- sort(unique(aliasmatrix[,2:3])))
    
    for(n in pnames){
        if(any(aliasmatrix[,2]==n)){
            p[n] <- aliasmatrix[aliasmatrix[,2]==n,1]
        }
        if(any(aliasmatrix[,3]==n)){
            p[n] <- p[n] + sum(1 - aliasmatrix[aliasmatrix[,3]==n,1])
        }
    }
    p<- p/length(p)
    p
}



ddiscrete <- function (x, probs, values = 1:length(probs))
{
    
    if (length(probs) != length(values))
        stop("ddiscrete: probs and values must have the same length.")
    if (sum(probs < 0) > 0)
        stop("ddiscrete: probs must not contain negative values.")
    if (!is.array(x) && !is.vector(x) && !is.factor(x))
        stop("ddiscrete: x must be an array or a vector or a factor.")
    
    p <- probs/sum(probs)
    
    y <- as.vector(x)
    l <- length(y)
    z <- rep(0,l)
    
    for (i in 1:l)
        if (any(values == y[i]))
            z[i] <- p[values == y[i]]
    
    z <- as.numeric(z)
    if (is.array(x))
        dim(z) <- dim(x)
    
    return(z)
}


pdiscrete <- function (q, probs, values = 1:length(probs))
{
    
    if (length(probs) != length(values))
        stop("pdiscrete: probs and values must have the same length.")
    if (sum(probs < 0) > 0)
        stop("pdiscrete: probs must not contain negative values.")
    if (!is.array(q) & !is.vector(q))
        stop("pdiscrete: q must be an array or a vector")
    
    p <- probs/sum(probs)
    
    y <- as.vector(q)
    l <- length(y)
    z <- rep(0,l)
    
    for (i in 1:l)
        z[i] <- sum(p[values <= y[i]])
    
    z <- as.numeric(z)
    if (is.array(q))
        dim(z) <- dim(q)
    
    return(z)
}

qdiscrete <- function (p, probs, values = 1:length(probs))
{
    
    if (length(probs) != length(values))
        stop("qdiscrete: probs and values must have the same length.")
    if (sum(probs < 0) > 0)
        stop("qdiscrete: probs must not contain negative values.")
    if (!is.array(p) & !is.vector(p))
        stop("qdiscrete: p must be an array or a vector")
    
    probs <- cumsum(probs)/sum(probs)
    
    y <- as.vector(p)
    l <- length(y)
    z <- rep(0,l)
    
    for (i in 1:l)
        z[i] <- length(values) - sum(y[i] <= probs) + 1
    
    z <- as.numeric(z)
    z <- values[z]
    if (is.array(p))
        dim(z) <- dim(p)
    
    return(z)
  }



