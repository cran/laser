### yule2rate
###
###  Likelihood of Branching Times Under Yule Model with 2 Speciation Rates
###
###  Author: Dan Rabosky <DLR32@cornell.edu>

yule2rate <- function(x, verbose = FALSE, ints = NULL, file = "out_yule2rate.txt") #3 parameters; 2 speciation rates
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
     #res <- res[res$LH != NA && res$LH != NaN]
     cat("----------------------------\nModel: yule 2 rate\n")
     cat("Log-likelihood and parameters:\n")
     cat("LH", summ$LH, "\nAIC", summ$aic, "\nr1", summ$r1, "\nr2", summ$r2)
     cat("\nst1", summ$st1, "\n")
     return(summ)
    
}
################################################################################

### yule3rate
###
### Likelihood of Branching Times Under Yule Model with 3 Speciation Rates

yule3rate <- function(x, ints = NULL, verbose = FALSE, file = "out_yule3rate.txt")
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
     #res <- res[res$LH != NA && res$LH != NaN]
     cat("----------------------------\nModel: yule 3 rate\n")
     cat("Log-likelihood and parameters:\n")
     cat("LH", summ$LH, "\nAIC", summ$aic, "\nr1", summ$r1, "\nr2", summ$r2, "\nr3", summ$r3)
     cat("\nst1", summ$st1, "\nst2", summ$st2, "\n")
     return(summ)
}
######################################################

### yule4rate
###
### Likelihood of Branching Times Under Yule Model with 4 Speciation Rates

yule4rate <- function(x, ints = NULL) #7 parameters: 4 speciation rates, 3 shift times
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
  
     cat("-----------------------------\n\n")
     cat("yule4rate:\nLog-likelihood & parameters:\n")
     cat("LH", res$LH, "\nAIC", res$aic, "\nr1", res$r1, "\nr2", res$r2, "\nr3", res$r3)
     cat("\nr4", res$r4)
     cat("\nst1", res$st1, "\nst2", res$st2, "\nst3", res$st3, "\n")
     
     
     return(res)


}
#########################################


### yule5rate
###
### Likelihood of Branching Times Under Yule Model with 5 Speciation Rates
###
### 


yule5rate <- function(x, ints = NULL)  #9 parameters with 5 speciation rates
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


     cat("-----------------------------\n\n")
     cat("yule5rate:\nLog-likelihood & parameters:\n")
     cat("LH", res$LH, "\nAIC", res$aic, "\nr1", res$r1, "\nr2", res$r2, "\nr3", res$r3)
     cat("\nr4", res$r4, "\nr5", res$r5)
     cat("\nst1", res$st1, "\nst2", res$st2, "\nst3", res$st3, "\nst4", res$st4, "\n")


     return(res)
     
}
########################################################