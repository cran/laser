
############ gamStat.batch ###################
###
###

gamStat.batch <- function(x, stat = NULL)
{

  if (!is.numeric(x)) stop("object x not of class 'numeric'")
  if (nrow(x) <= 1)
  {
    print("Error in input data format\n")
    stop()
  }
  res <- list()
  temp <- list()

  for (i in 1:nrow(x))
  {
    res$gamstat[i] <- batchgam(x[i,])
  }

  v <- res$gamstat
  v <- as.vector(v)
  v <- sort(v)

  gamcrit <- v[.05*length(v)]
  cat("------------------------\n")
  cat("Critical value of gamma statistic (alpha = 0.05):", gamcrit, "\n")
  if (!is.null(stat))
  {
    tvec <- v[v < stat]
    pval <- length(tvec)/length(v)
    cat("\nProb(x | Ho: rates have not decreased):", pval, "\n\n")
  }
  res$gamcrit <- gamcrit
  if (!is.null(stat))
    res$pval <- pval
  return(res)

 }

######        'private' functions       #######
###
### batchgam
###
### Used to Calculate Gamma Statistic for Batch of Phylogenies
###
###

batchgam <- function(x)
{

    res <- list()
    x <- rev(sort(x))
    N <- length(x)+1
    b <- sort(x)
    z <- rev(c(b[1], diff(b)))
    T <- sum((2:N)*z)
    s1 <- (1/(N-2) * sum(cumsum((2:(N-1))*z[1:(N-2)])) - T/2)
    res = s1/(T*sqrt(1/(12*(N-2))))

    return(res)
}
