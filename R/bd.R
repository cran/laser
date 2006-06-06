### bd
###
###  Fit Constant Rate Birth-Death Model to Set of Branching Times
###

bd <- function(x, ai = c(0.1, 0.5, 0.9)) #new 'optim' version, 6.4.06
{
  N <- length(x)+1
  b <- sort(x)
  z <- rev(c(b[1], diff(b)))
  x <- c(0, x)
  res <- list()

  mlbd <- function(v)
  {
    r <- v[1]
    a <- v[2]
    -( sum(log(1:(N-1))) + ((N-2)*log(r))
        + (r*sum(x[3:N]))
        +(N*log(1-a)) - 2 * sum(log(exp(r * x[2:N])-a)))
  }

  for (k in 1:length(ai))
  {
    temp <- suppressWarnings(optim(c(.2, ai[k]), mlbd))
    if (temp$par[2] <= 0)
    {
      temp <- IpureBirth(x[2:length(x)])
      if (k == 1 || (k > 1 && temp$LH > res$LH))
      {
        res$LH <- temp$LH
        res$r1 <- temp$r1
        res$a <- 0
      }
    }
    else if (k == 1 || (k > 1 && res$LH < -temp$value))
    {
        res$LH <- -temp$value
        res$r1 <- temp$par[1]
        res$a <- temp$par[2]
    }
  }
  res$aic <- (-2*res$LH) + 4
  cat("----------------------------\n")
  cat("Model: bd (rate-constant)\n\n")
  cat("Log-likelihood and parameters:\nLH:", res$LH, "\n")
  cat("AIC:", res$aic, "\nr = (S - E) =", res$r1, "\n")
  cat("a = (E/S) =", res$a, "\n\n")
  return(res)

}
