lca <- function(x, k, niter=100, matchdata=FALSE, verbose=FALSE)
{

    ## if x is a data matrix -> create patterns
    if (is.matrix(x))
    {
        if (matchdata)
        {
            x <- countpattern(x, matching=TRUE)
            xmat <- x$matching
            x <- x$pat
        }
        else
            x <- countpattern(x, matching=FALSE)
    }
    else   ## if no data ist given, matchdata must be FALSE
        matchdata <- FALSE

    n <- sum(x)
    npat <- length(x)
    nvar <- round(log(npat)/log(2))
    
    ## build matrix of all possible binary vectors
    b <- matrix(0, 2^nvar, nvar)
    for (i in 1:nvar)
        b[, nvar+1-i] <- rep(rep(c(0,1),c(2^(i-1),2^(i-1))),2^(nvar-i))

    ## initialize probabilities
    classprob <- runif(k)
    classprob <- classprob/sum(classprob)
    names(classprob) <- 1:k
    p <- matrix(runif(nvar*k), k)

    
    pas <- matrix(0, k, npat)
    classsize <- numeric(k)
    
    for (i in 1:niter)
    {
        for (j in 1:k)
        {
            ## P(pattern|class)
            mp <- t(b)*p[j,]+(1-t(b))*(1-p[j,])
            pas[j,] <- drop(exp(rep(1,nvar)%*%log(mp))) # column product
        }
        ##  P(pattern|class)*P(class)
        pas <- t(t(pas)*classprob)        
        
        ## P(class|pattern)
        sump <- drop(rep(1,k)%*%pas)  # column sums
        pas <- t(t(pas)/sump)

        spas <- t(t(pas)*x)
        classsize <- drop(spas%*%rep(1,npat))  # row sums
        classprob <- classsize/n
        p <- pas%*%(x*b)/classsize
        if (verbose)
            cat("Iteration:", i, "\n")
    }

    for (j in 1:k)
    {
        mp <- t(b)*p[j,]+(1-t(b))*(1-p[j,])
        pas[j,] <- drop(exp(rep(1,nvar)%*%log(mp)))*classprob[j]
                                        # column product
    }

    ## LogLikelihood
    pmust <-  drop(rep(1,k)%*%pas)  # column sums
    ll <- sum(x*log(pmust))

    ## Likelihoodquotient
    xg0 <- x[x>0]
    ll0 <- sum(xg0*log(xg0/n))
    lq <- 2*(ll0-ll)

    ## bic
    bic <- -2*ll+log(n)*(k*(nvar+1)-1)
    bicsat <- -2*ll0+log(n)*(2^nvar-1)
    
    ## chisq
    ch <- sum((x-n*pmust)^2/(n*pmust))
    
    ## P(class|pattern)
    sump <- drop(rep(1,k)%*%pas)  # column sums
    pas <- t(t(pas)/sump)

    mat <- max.col(t(pas))
    if (matchdata)
        mat <- mat[xmat]

    colnames(p) <- 1:nvar
    rownames(p) <- 1:k
    
    lcaresult <- list(classprob=classprob, p=p, matching=mat,
                      logl=ll, loglsat=ll0,
                      chisq=ch, lhquot=lq, bic=bic, bicsat=bicsat, n=n,
                      np=(k*(nvar+1)-1), matchdata=matchdata)

    class(lcaresult) <- "lca"
    return(lcaresult)
}


print.lca <- function(l)
{
    cat("LCA-Result\n")
    cat("----------\n\n")

    cat("Datapoints:", l$n, "\n")
    cat("Classes:   ", length(l$classprob), "\n")
    cat("Probability of classes\n")
    print(round(l$classprob,3))

    cat("Itemprobabilities\n")
    print(round(l$p,2))
}

summary.lca <- function(l)
{
    nvar <- ncol(l$p)
    l$npsat <- 2^nvar-1
    l$df <- 2^nvar-1-l$np
    l$pvallhquot <- 1-pchisq(l$lhquot,l$df)
    l$pvalchisq <- 1-pchisq(l$chisq,l$df)
    l$k <- length(l$classprob)

    ## remove unnecessary list elements
    l$classprob <- NULL
    l$p <- NULL
    l$matching <- NULL
    
    class(l) <- "summary.lca"
    return(l)
}


print.summary.lca <- function(l)
{
    cat("LCA-Result\n")
    cat("----------\n\n")

    cat("Datapoints:", l$n, "\n")
    cat("Classes:   ", l$k, "\n")

    cat("\nGoodness of fit statistics:\n\n")
    cat("Number of parameters, estimated model:", l$np, "\n")
    cat("Number of parameters, saturated model:", l$npsat, "\n")
    cat("Log-Likelihood, estimated model:      ", l$logl, "\n")
    cat("Log-Likelihood, saturated model:      ", l$loglsat, "\n")

    cat("\nInformation Criteria:\n\n")
    cat("BIC, estimated model:", l$bic, "\n")
    cat("BIC, saturated model:", l$bicsat, "\n")

    cat("\nTestStatistics:\n\n")
    cat("Likelihood ratio:  ", l$lhquot,
        "  p-val:", l$pvallhquot, "\n")
    cat("Pearson Chi^2:     ", l$chisq,
        "  p-val:", l$pvalchisq, "\n")
    cat("Degress of freedom:", l$df, "\n")
}


bootstrap.lca <- function(l, nsamples=10, lcaiter=30, verbose=FALSE)
{

    n <- l$n
    classprob <- l$classprob
    nclass <- length(l$classprob)
    prob <- l$p
    nvar <- ncol(l$p)
    npat <- 2^nvar

    ## build matrix of all possible binary vectors
    b <- matrix(0, npat, nvar)
    for (i in 1:nvar)
        b[, nvar+1-i] <- rep(rep(c(0,1),c(2^(i-1),2^(i-1))),2^(nvar-i))

    ll <- lq <- ll0 <- ch <- numeric(nsamples)
    
    for (i in 1:nsamples)
    {
        ## generate data
        cm <- sample(1:nclass, size=n, replace=TRUE, prob=classprob)
        x <- matrix(runif(n*nvar), nrow=n)
        x <- (x<prob[cm,])+0
        x <- countpattern(x, matching=FALSE)

        if (verbose)
            cat ("Start of Bootstrap Sample", i, "\n")
        lc <- lca(x, nclass, niter=lcaiter, verbose=verbose)
        
        ll[i] <- lc$logl
        ll0[i] <- lc$loglsat
        lq[i] <- lc$lhquot
        ch[i] <- lc$chisq

        if (verbose)
            cat("LL:", ll[i], " LLsat:", ll0[i],
                " Ratio:", lq[i], " Chisq:", ch[i], "\n")
    }


    lratm <- mean(lq)
    lrats <- sd(lq)
    chisqm <- mean(ch)
    chisqs <- sd(ch)

    zl <- (l$lhquot-lratm)/lrats
    zc <- (l$ch-chisqm)/chisqs
    pzl <- 1-pnorm(zl)
    pzc <- 1-pnorm(zc)
    pl <- sum(l$lhquot<lq)/length(lq)
    pc <- sum(l$ch<ch)/length(ch)
    
    lcaboot <- list(logl=ll, loglsat=ll0, lratio=lq, 
                    lratiomean=lratm, lratiosd=lrats,
                    lratioorg=l$lhquot, zratio=zl,
                    pvalzratio=pzl, pvalratio=pl,
                    chisq=ch, chisqmean=chisqm, chisqsd=chisqs,
                    chisqorg=l$ch, zchisq=zc,
                    pvalzchisq=pzc, pvalchisq=pc,
                    nsamples=nsamples, lcaiter=lcaiter)
    class(lcaboot) <- "bootstrap.lca"
    return(lcaboot)
}

print.bootstrap.lca <- function(bl)
{
    cat("Bootstrap of LCA\n")
    cat("----------------\n\n")

    cat ("Number of Bootstrap Samples:    ", bl$nsamples, "\n")
    cat ("Number of LCA Iterations/Sample:", bl$lcaiter, "\n")
    
    cat("Likelihood Ratio\n\n")
    cat("Mean:", bl$lratiomean, "\n")
    cat("SDev:", bl$lratiosd, "\n")
    cat("Value in Data Set:", bl$lratioorg, "\n")
    cat("Z-Statistics:     ", bl$zratio, "\n")
    cat("P(Z>X):           ", bl$pvalzratio, "\n")
    cat("P-Val:            ", bl$pvalratio, "\n\n")

    cat("Pearson's Chisquare\n\n")
    cat("Mean:", bl$chisqmean, "\n")
    cat("SDev:", bl$chisqsd, "\n")
    cat("Value in Data Set:", bl$chisqorg, "\n")
    cat("Z-Statistics:     ", bl$zchisq, "\n")
    cat("P(Z>X):           ", bl$pvalzchisq, "\n")
    cat("P-Val:            ", bl$pvalchisq, "\n\n")
}

predict.lca <- function(lcares, x)
{
    if (lcares$matchdata)
        stop("predict.lca: only possible, if lca has been called with matchdata=FALSE")
    else
    {
        x <- countpattern(x, matching=TRUE)
        return(lcares$matching[x$matching])
    }
}
