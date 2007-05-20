`yule4rate` <-
function(x, ints = NULL) #7 parameters: 4 speciation rates, 3 shift times
{
  #uses observed shift times only; computational time very high otherwise#
  if (!is.numeric(x)) stop("object x not of class 'numeric'")
  N <- length(x)+1
  if (is.null(ints))
  {
    tv <- x[2:(N-2)]
    tv <- unique(tv)
    stvec <- combinations(length(tv), 3, tv)
  }
  else
  {
    inc <- (x[2] - x[length(x)])/ints
    iv <- seq(x[2], length.out = ints, by = -inc)
    stvec <- combinations(length(iv), 3, iv)
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
     v4 <- 0
      v1 <- yuleint2(x, x[1], stvec[i, 1])
      v2 <- yuleint2(x, stvec[i, 1], stvec[i, 2])
      v3 <- yuleint2(x, stvec[i,2], stvec[i, 3])
      v4 <- yuleint2(x, stvec[i,3], 0)
      LHtemp <- v1$LH + v2$LH + v3$LH + v4$LH
      
      if (is.finite(LHtemp) == TRUE && (i == 1 || (i > 1 && LHtemp > res$LH)))
      {
      res$LH <- v1$LH + v2$LH + v3$LH + v4$LH
      res$st1 <- stvec[i, 1]
      res$st2 <- stvec[i, 2]
      res$st3 <- stvec[i, 3]
      res$r1 <- v1$smax
      res$r2 <- v2$smax
      res$r3 <- v3$smax
      res$r4 <- v4$smax
      }
      
  }
     res$aic <- (-2*res$LH)+ 14
     res <- as.data.frame(res)
     res <- na.omit(res)

     
     return(res)


}

