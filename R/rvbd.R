### rvbd
###
### Likelihood of Branching Times Under Birth-Death Model with Rate Shift
###

rvbd <- function(x, ai = c(0.1, .5, .95), ints = NULL, verbose = FALSE, file = "out_rvbd.txt")
{
    checkbasal(x)
    if (!is.numeric(x)) stop("object x not of class 'numeric'")
    x <- rev(sort(x))
    res <- data.frame()
    temp <- data.frame()
    N <- length(x)+1
    Nvec <- 2:N
    x <- c(0, x)
    rlist <- list()
    if (is.null(ints))
      stvec <- x[4:(length(x)-2)]
    else
      stvec <- seq(x[4], x[length(x)-2], length.out = ints)
    LHra  <- function(v)
    {
    r1 <- v[1]
    r2 <- v[2]
    a <- v[3]

   -(sum(log(1:(N-1))) + (i-2)*log(r1) + (N - i)*log(r2)
    + sum((x[3:i]-st)*r1 + st*r2) + sum(x[(i+1):N]*r2)
    + N*log(1-a)
    - 2 * sum(log((exp((x[2:i]-st)*r1 + (st*r2)) - a)))
    - 2 * sum(log((exp(x[(i+1):N]*r2) - a))))
    }
    for (j in 1:length(stvec))
    {
      for (z in 2:(length(x)-1))
      {
        if (x[z] >= stvec[j] && x[z+1] < stvec[j])
        {
          i = z
        }
      }

      #go through different initial parameterizations for 'a'
      for (k in 1:length(ai))
      {
        st <- stvec[j]
        temp <- suppressWarnings(optim(c(.3, .3, ai[k]), LHra))

        if (temp$par[3] < 0)
        {
          t1 <- yuleint2(x[2:length(x)], x[2], st)
          t2 <- yuleint2(x[2:length(x)], st, 0)
          temp$par[1] <- t1$smax
          temp$par[2] <- t2$smax
          temp$par[3] <- 0
          temp$value <- -(t1$LH + t2$LH)
        }
        if (k == 1)
        {
          res <- temp
        }
        else if ((k > 1) && (temp$value < res$value))
        {
          res <- temp
        }

      }
    rlist$r1[j] <- res$par[1]
    rlist$r2[j] <- res$par[2]
    rlist$LH[j] <- -res$value
    rlist$st[j] <- st
    rlist$a[j] <- res$par[3]
    }
    rlist <- as.data.frame(rlist)
    rlist <- na.omit(rlist)
    if (verbose == TRUE)
    {
      sink(file)
      cat("i  LH  r1  r2  a st\n")
      for (i in 1:length(stvec))
        cat(i, rlist$LH[i], rlist$r1[i], rlist$r2[i], rlist$a[i], rlist$st[i], "\n")
      sink()
    }

    summ <- rlist[rlist$LH == max(rlist$LH), ]
    summ$aic <- (-2*summ$LH) + 8
    cat("------------------------\nModel: Rate Variable Birth-Death: r1, r2, a, ts\n")
    cat("Log-likelihood and parameter estimates\n")
    cat("LH", summ$LH, "\nAIC", summ$aic, "\nr1", summ$r1, "\nr2", summ$r2, "\na", summ$a, "\nst", summ$st, "\n")

    return(summ)
}

#######################################################