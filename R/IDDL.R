`IDDL` <-
function(x)
{
  #calculates likelihoods under DD logistic model
  N <- length(x)+1
  b <- sort(x)
  z <- rev(c(b[1], diff(b)))
  res <-list()
  #negative log-LH function to be minimized
  ddfunc <- function(r, k)
  {
    -(sum(log(2:(N-1))) + (N-2)*log(r) + sum(log(1-((2:(N-1))/k)))
    - sum((2:N)*r*z) + sum(z*r*(2:N)^2)/k)
  }
   temp <- nlm(function(p) ddfunc(p[1], p[2]), c(.5, N*2), hessian = TRUE)
   res$LH <- -temp$minimum
   res$r1 <- temp$estimate[1]
   res$k <- temp$estimate[2]
   return(res)
}

