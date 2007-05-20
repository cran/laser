`yule5rate` <-
function(x, ints = NULL)  #9 parameters with 5 speciation rates
{
  #uses observed shift times only#
  #run times can be high...
  if (!is.numeric(x)) stop("object x not of class 'numeric'")
  N <- length(x)+1
  
  if (is.null(ints))
  {
    tv <- x[2:(N-2)]
    tv <- unique(tv)
    stvec <- combinations(length(tv), 4, tv)
  }
  else
  {
    inc <- (x[2] - x[length(x)])/ints
    iv <- seq(x[2], length.out = ints, by = -inc)
    stvec <- combinations(length(iv), 4, iv)
   }
  for (i in 1:(length(stvec[,1])))
  {
    stvec[i,] <- rev(sort(stvec[i,]))
  }
  
  
  
  res <- list()
  #if (verbose == TRUE)
  #{
   # sink(file)
    #cat("i  LH  r1  r2  r3  r4  r5  st1 st2 st3 st4\n")
  #}
  for (i in 1:(length(stvec[,1])))
  {

     v1 <- 0
     v2 <- 0
     v3 <- 0
     v4 <- 0
     v5 <- 0
      v1 <- yuleint2(x, x[1], stvec[i, 1])
      v2 <- yuleint2(x, stvec[i, 1], stvec[i, 2])
      v3 <- yuleint2(x, stvec[i,2], stvec[i, 3])
      v4 <- yuleint2(x, stvec[i,3], stvec[i ,4])
      v5 <- yuleint2(x, stvec[i, 4], 0)
      LHtemp <- v1$LH +v2$LH + v3$LH + v4$LH +v5$LH
      if (is.finite(LHtemp) == TRUE && (i == 1 || (i > 1 && LHtemp > res$LH)))
      {
      res$r1 <- v1$smax
      res$r2 <- v2$smax
      res$r3 <- v3$smax
      res$r4 <- v4$smax
      res$r5 <- v5$smax
      res$LH <- LHtemp
      res$st1 <- stvec[i, 1]
      res$st2 <- stvec[i, 2]
      res$st3 <- stvec[i, 3]
      res$st4 <- stvec[i, 4]
      }
      #if (verbose == TRUE)
      #  cat(i, res$LH[i], res$r1[i], res$r2[i], res$r3[i], res$r4[i], res$r5[i], res$st1[i],
      #    res$st2[i], res$st3[i], res$st4[i], "\n")
  }
     #if (verbose == TRUE)
     # sink()
     res$aic <- (-2 * res$LH) + 18
     res <- as.data.frame(res)
     res <- na.omit(res)

     return(res)
     
}

