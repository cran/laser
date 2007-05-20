`yule2rate` <-
function(x, verbose = FALSE, ints = NULL, file = "out_yule2rate.txt") #3 parameters; 2 speciation rates
{
    checkbasal(x)
    if (!is.numeric(x)) stop("object x not of class 'numeric'")
    x <- rev(sort(x))
    if (is.null(ints))
      stvec <- x[3:(length(x)-2)]
    else
      stvec <- seq(x[3], x[length(x)-2], length.out = ints)
    #stvec is a vector of possible shift times, determined either by
    # observed branching times or by specified intervals.

    res <- list()
    #verbose TRUE prints LH, r1, r2 for each shift rate.
    if (verbose == TRUE)
      cat("i",  "LH",  "r1",  "r2",  "st1\n", file = file, sep = "\t")
    
    for (i in 1:(length(stvec)))
    {
      v1 <- yuleint2(x, x[1], stvec[i])
      v2 <- yuleint2(x, stvec[i], 0)
      res$LH[i] <- v1$LH + v2$LH
      res$st1[i] <- stvec[i]
      res$r1[i] <- v1$smax
      res$r2[i] <- v2$smax
      if (verbose == TRUE)
        cat(i, res$LH[i], res$r1[i], res$r2[i], res$st1[i], "\n", file = file, append = TRUE, sep = "\t")
    }
     res <- as.data.frame(res)
     res <- na.omit(res)
     summ <- res[res$LH == max(res$LH), ]
     summ$aic <- (-2*summ$LH) + 6

     return(summ)
    
}

