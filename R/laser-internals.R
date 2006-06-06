
### internal functions for laser ##################
###
### All functions written by Dan Rabosky <DLR32@cornell.edu>
### except for 'combinations', which was modified by
### Gregory R. Warnes <gregory.r.warnes@pfizer.com> and made available 
### as part of the 'gtools' package.  Original is from email 
### by Brian D Ripley <ripley@stats.ox.ac.uk> to r-help
### dated Tue, 14 Dec 1999 11:14:04 +0000 (GMT) in response to
### Alex Ahgarin <datamanagement@email.com>.  Original version 
### of 'combinations' was named "subsets" and was Written by Bill Venables.


### Ibd
###
### bd modelfitting for fitdAICrc

Ibd <- function(x) #new 'optim' version, 6.4
{
  N <- length(x)+1
  b <- sort(x)
  z <- rev(c(b[1], diff(b)))
  x <- c(0, x)
  res <- list()
  ai <- c(.1, .9)
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

  return(res)

}
### Irvbd
###
### rvbd model-fitting for fitdAICrc
Irvbd <- function(x, ai = c(.1, .5, .9), ints = NULL)
{
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

   LHra  <- function(v)   #v[1], v[2], v[3] v[4]
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
        res <- temp
      else if ((k > 1) && (temp$value < res$value))
        res <- temp
    }
   rlist$r1[j] <- res$par[1]
   rlist$r2[j] <- res$par[2]
   rlist$LH[j] <- -res$value
   rlist$st[j] <- st
   rlist$a[j] <- res$par[3]
  }
  res2 <- list()
  res2$LH <- max(rlist$LH)
  res2$r1 <- rlist$r1[rlist$LH == max(rlist$LH)]
  res2$r2 <- rlist$r2[rlist$LH == max(rlist$LH)]
  res2$a <- rlist$a[rlist$LH == max(rlist$LH)]
  res2$st <- rlist$st[rlist$LH == max(rlist$LH)]
  return(res2)

}
### IDDX
###
### DDX model-fitting for fitdAICrc

IDDX <- function(x)
{

  N <- length(x)+1
  b <- sort(x)
  z <- rev(c(b[1], diff(b)))
  #initial parameterization for nlm fx: estimates under yule model
  s1 <- IpureBirth(x)
  #print(s1$r1)
  res <- list()
  temp <- list()
  ddfunc <- function(r, v)
  {
    -(sum(log(2:(N-1))) + (N-2)*log(r) - sum(v*log(2:(N-1)))
    - sum(((2:N)^(1-v))*z *r) )
  }
  temp <- nlm(function(p) ddfunc(p[1], p[2]), c(s1$r1, 0), hessian = TRUE)

  res$LH <- -temp$minimum
  res$r1 <- temp$estimate[1]
  res$xp <- temp$estimate[2]

   return(res)
}

### IDDL
###
### DDL model-fitting for fitdAICrc

IDDL <- function(x)
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

### IpureBirth
###
### pureBirth modelfitting for fitdAICrc

IpureBirth <- function(x)
{
  res <- list()
  temp <- yuleint2(x, x[1], 0)
  res$LH <- temp$LH
  res$r1 <- temp$smax
  return(res)
}

### Iyule2rate
###
### yule2rate modelfitting function for fitdAICrc

Iyule2rate <- function(x, ints = NULL)
{
 res <- list()
 temp <- list()

 if (is.null(ints))
      stvec <- x[3:(length(x)-2)]
 else
      stvec <- seq(x[3], x[length(x)-2], length.out = ints)

 for (i in 1:length(stvec))
 {
    v1 <- yuleint2(x, x[1], stvec[i])
    v2 <- yuleint2(x, stvec[i], 0)
    temp$LH[i] <- v1$LH + v2$LH
    temp$st[i] <- stvec[i]
    temp$r1[i] <- v1$smax
    temp$r2[i] <- v2$smax
 }
    res$LH <- max(temp$LH)
    res$st <- temp$st[temp$LH == max(temp$LH)]
    res$r1 <- temp$r1[temp$LH == max(temp$LH)]
    res$r2 <- temp$r2[temp$LH == max(temp$LH)]
    return (res)
}

### Iyule3rate
###
### yule3rate modelfitting for fitdAICrc

Iyule3rate <- function(x, ints = NULL)
{
  x <- rev(sort(x))
  N <- length(x)+1

  if (is.null(ints))
  {
    tv <- x[2:(N-2)]
    tv <- unique(tv)
    stvec <- combinations(length(tv), 2, tv)
  }
  else
  {
    inc <- (x[2] - x[length(x)])/ints
    iv <- seq(x[2], length.out = ints, by = -inc)
    stvec <- combinations(length(iv), 2, iv)
   }
  for (i in 1:(length(stvec[,1])))
  {
    stvec[i,] <- rev(sort(stvec[i,]))
  }

  res <- list()
  for (i in 1:(length(stvec[,1])))
  {
     v1 <- 0
     v2 <- 0
     v3 <- 0
      v1 <- yuleint2(x, x[1], stvec[i, 1])
      v2 <- yuleint2(x, stvec[i, 1], stvec[i, 2])
      v3 <- yuleint2(x, stvec[i,2], 0)
      res$LH[i] <- v1$LH + v2$LH + v3$LH
      res$r1[i] <- v1$smax
      res$r2[i] <- v2$smax
      res$r3[i] <- v3$smax
      res$st[i] <- stvec[i, 1]
      res$st2[i] <- stvec[i, 2]

  }

     res <- as.data.frame(res)
     res <- na.omit(res)
     summ <- res[res$LH == max(res$LH), ]
     return(summ)
}



### combinations
###
### Combinations of size n from v
###
### From email by Brian D Ripley <ripley@stats.ox.ac.uk> to r-help
### dated Tue, 14 Dec 1999 11:14:04 +0000 (GMT) in response to
### Alex Ahgarin <datamanagement@email.com>.  Original version was
### named "subsets" and was Written by Bill Venables.
###

combinations <- function(n, r, v = 1:n, set = TRUE, repeats.allowed=FALSE)
{
  if(mode(n) != "numeric" || length(n) != 1
     || n < 1 || (n %% 1) != 0) stop("bad value of n")
  if(mode(r) != "numeric" || length(r) != 1
     || r < 1 || (r %% 1) != 0) stop("bad value of r")
  if(!is.atomic(v) || length(v) < n)
    stop("v is either non-atomic or too short")
  if( (r > n) & repeats.allowed==FALSE)
    stop("r > n and repeats.allowed=FALSE")
  if(set) {
    v <- unique(sort(v))
    if (length(v) < n) stop("too few different elements")
  }
  v0 <- vector(mode(v), 0)
  ## Inner workhorse
  if(repeats.allowed)
    sub <- function(n, r, v)
      {
        if(r == 0) v0 else
        if(r == 1) matrix(v, n, 1) else
        if(n == 1) matrix(v, 1, r) else
        rbind( cbind(v[1], Recall(n, r-1, v)),
              Recall(n-1, r, v[-1]))
      }
  else
    sub <- function(n, r, v)
      {
        if(r == 0) v0 else
        if(r == 1) matrix(v, n, 1) else
        if(r == n) matrix(v, 1, n) else
        rbind(cbind(v[1], Recall(n-1, r-1, v[-1])),
              Recall(n-1, r, v[-1]))
      }
  sub(n, r, v[1:n])
}

###  yuleint2
###
### Internal workhorse for various functions
###

yuleint2 <- function(x, st1, st2)
{
  nvec <- 2:(length(x)+1)
  nv <- x[(x < st1) & (x >= st2)]
  lo <- max(nvec[x >= st1])
  up <- max(nvec[x >= st2])
  
  res <- list()
  if (st1 <= x[1])
    {nv <- c(st1, nv) - st2}
  else
    {nv <- nv - st2}

  smax <- (up-lo)/(lo*nv[1] + sum(nv[2:(up-lo +1)]))

  s1 <- sum(log(lo:(up-1)))
  s2 <- (up-lo)*log(smax)
  s3 <- -(lo*nv[1] + sum(nv[2:(up-lo+1)]))*smax
  res$smax <- smax
  res$LH <- s1 + s2 + s3

  return(res)

}

###############################################

checkbasal <- function(x)
{
  if (x[1] == x[2])
    stop("Invalid data format: basal polytomy")
}