`getLambda.internal` <-
function (zmat, rootnode, rbounds, para = 0.01, eps, combined = TRUE) 
{
    int <- zmat[zmat[, 2] > rootnode, ]
    term <- zmat[zmat[, 2] < rootnode, ]
    nint <- nrow(int)
    nterm <- nrow(term)
    betaF <- function(r, t1) {
        xf <- (exp(r * t1) - 1)/(exp(r * t1) - eps)
        xf
    }
    Lfunc_tax <- function(p) {
        r <- p
        (sum(log(1 - betaF(r, term[1:nterm, 4]))) + sum((term[1:nterm, 
            5] - 1) * log(betaF(r, term[1:nterm, 4]))))
    }
    Lfunc_phy <- function(p) {
        r <- p
        (nint * log(r) - r * sum(int[1:nint, 4]) - sum(log(1 - 
            (eps * exp(-r * int[1:nint, 3])))))
    }
    Lfunc_comb <- function(p) {
        r <- p
        (sum(log(1 - betaF(r, term[1:nterm, 4]))) + sum((term[1:nterm, 
            5] - 1) * log(betaF(r, term[1:nterm, 4]))) + nint * 
            log(r) - r * sum(int[1:nint, 4]) - sum(log(1 - (eps * 
            exp(-r * int[1:nint, 3])))))
    }

    res <- list()
    if (combined == TRUE) {
        if (nrow(int) == 0) 
            tempres <- optimize(Lfunc_tax, interval = rbounds, maximum = TRUE)
        else tempres <- optimize(Lfunc_comb, interval = rbounds, maximum = TRUE)
    }
    else {
        tempres <- optimize(Lfunc_tax, interval = rbounds, 
            maximum = TRUE)
    }
    res$LH <- tempres$objective
    res$lambda <- tempres$maximum/(1 - eps)
    res$r <- tempres$maximum
    res$eps <- eps
    res <- as.data.frame(res)
    return(res)
}

