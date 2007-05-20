`yule3rate` <-
function(x, ints = NULL, verbose = FALSE, file = "out_yule3rate.txt")
{
  if (!is.numeric(x)) stop("object x not of class 'numeric'")
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
  if (verbose == TRUE)
    cat("i",  "LH",  "r1",  "r2", "r3", "st1", "st2\n", file = file, sep = "\t")
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
      res$st1[i] <- stvec[i, 1]
      res$st2[i] <- stvec[i, 2]
      # LH estimates for each temporal window - uncomment to output...
      #res$L1[i] <- v1$LH
      #res$L2[i] <- v2$LH
      #res$L3[i] <- v3$LH
      if (verbose == TRUE)
        cat(i, res$LH[i], res$r1[i], res$r2[i], res$r3[i], res$st1[i], res$st2[i], "\n", file = file, append = TRUE, sep = "\t")
  }
     res <- as.data.frame(res)
     res <- na.omit(res)
     summ <- res[res$LH == max(res$LH), ]
     summ$aic <- (-2*summ$LH) + 10

     return(summ)
}

