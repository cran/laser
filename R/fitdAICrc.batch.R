

fitdAICrc.batch <- function(x, modelset = c("pureBirth", "bd", "DDL", "DDX", "yule2rate"),
                verbose = TRUE, file = "out_daic.batch.txt", ints = NULL,
                alternative = NULL, stat = NULL)
{
  if (!is.numeric(x)) stop("object x not of class 'numeric'")
  if (nrow(x) <= 1)
  {
    print("Error in input data format\n")
    stop()
  }
  res <- list()
  temp <- list()

  #batch dAICrc statistic
  setvec <- c("pureBirth", "bd", "DDX", "DDL", "rvbd", "yule2rate", "yule3rate") 
    #vector of all possible candidate models
  abbvec <- c("pb", "bd", "DX", "DL", "rvbd", "y2r8", "y3r8")
    #abbreviations of models
  
  # if verbose TRUE, prints header of variable names for output file
  if (verbose == TRUE)
  {
    sink(file)
    cat("bestmodel\tdAICrc\tsp\t")
    pvec <- c("LH", "AIC", "r1", "r2", "r3", "a", "x", "k", "st", "st2")
    for (i in 1:length(modelset))
    {
      tname <- paste(abbvec[setvec == modelset[i]], "_", pvec, "\t", sep = "")
      cat(tname)
    }
    cat("\n")
    sink()
  }

  #sends each set of branching times to dAICrc.batch, which returns summary.
  
  for (i in 1:nrow(x))
  {
    if ((i %% 10) == 0)
      cat("set", i, "\n")
    xvec <- as.numeric(x[i,])
    temp <- dAICrc.batch(xvec, modelset = modelset, verbose = verbose, file = file, ints = ints)
    res <- rbind(res, temp)
  }
  
  # summarizing results:
  if (verbose == TRUE)
    summarydAICrc(file = file, alternative = alternative, stat = stat)

  return(res)
}



### summarydAICrc
###
###

summarydAICrc <- function(file = NULL, alternative = NULL, stat = NULL)
{
  # 'file' is output file fitdAICrc.batch
  
  # stat is calculated dAICrc from real data that you want to compare to null distribution

  if (is.null(file))
    stop("Cannot read file for 'verbose' output from the batch run\n")
  x <- read.table(file = file, header = TRUE)
  x <- t(na.omit(t(x)))
  x <- as.data.frame(x)
  #preceding lines take output file from batchPhyStat and delete all 'NA' columns

  vd <- x$dAICrc[x$sp == "RD"]  #vector of all dAICrc's where rates DECREASED over time
  vin <- x$dAICrc[x$sp == "RI"]
    #vector of all dAICrc's where rates INCREASED over time
  vnc <- x$dAICrc[x$sp == "NC"]  #all dAICrc's with RATE CONSTANT as best model
  vd <- as.vector(vd)
  vin <- as.vector(vin)
  vnc <- as.vector(vnc)
  vd <- as.numeric(vd)
  vin <- as.numeric(vin)
  vnc <- as.numeric(vnc)
  #not sure why I have to do this in so many steps, but otherwise does not work

  ### one-tailed test for rate DECREASE
  cat("------------------------------------------\n")
  if (!is.null(stat))
    cat("dAICrc calculated: ", stat, "\n")
  cat("Summary:\n")
  if (!is.null(alternative) && alternative == "decrease")
  {
    vd <- sort(vd)
    v <- c(vnc, vin, vd)

    dAICcrit_1TD <- v[.95 * length(v)]
    if (!is.null(stat))
    {
      temp <- vd[vd > stat]
      pd <- length(temp)/length(v)
      if (length(temp) == length(vd))
        sgn <- ">"
      else if (length(temp) == 0)
        sgn <- "="
      else
        sgn <- "<"
    }
    cat("Test for rate decrease over time:\n")
    cat("dAICrc_crit (one-tailed, alpha = 0.05):  ", dAICcrit_1TD, "\n");
    if (!is.null(stat))
      cat("prob under Ho (no decrease in rates: ", paste(sgn), pd, "\n\n")
  }
  #one tailed test for rate increase
  if (!is.null(alternative) && alternative == "increase")
  {
    vin <- sort(vin)
    v <- c(vnc, vd, vin)
    dAICcrit_1TI <- v[.95 * length(v)]
    if (!is.null(stat))
    {
      temp <- vin[vin > stat]
      pd <- length(temp)/length(v)
      if (length(temp) == length(vin))
        sgn <- ">"
      else if (length(temp) == 0)
        sgn <- "="
      else
        sgn <- "<"
    }
    cat("Test for rate increase over time:\n")
    cat("dAICrc_crit (one-tailed, alpha = 0.05):  ", dAICcrit_1TI, "\n");
    if (!is.null(stat))
      cat("prob under Ho (no increase in rates: ", paste(sgn), pd, "\n\n")
  }
  #two tailed test for rate variation
  if (is.null(alternative))
  {
    v <- c(vnc, vd, vin)
    v <- sort(v)
    dAICcrit_2T <- v[.95 * length(v)]
    if (!is.null(stat))
    {
      temp <- v[v > stat]
      pd <- length(temp)/length(v)
      print(pd)
      if (length(temp) == length(v))
        sgn <- ">"
      else if (length(temp) == 0)
        sgn <- "="
      else
        sgn <- "<"
    }
    cat("Test for rate constancy:\n")
    cat("dAICrc crit (2-tailed, alpha = 0.05): ", dAICcrit_2T, "\n")
    if (!is.null(stat))
      cat("prob (x | rate constant): ", paste(sgn), pd, "\n")
  }
  return()
}

### dAICrc.batch
###
### Required function for generating null distribution of dAICrc

dAICrc.batch <- function(x, modelset = mlist, verbose2 = TRUE, file, ints)
{
  res <- list()
  temp <- list()

  for (i in 1:length(modelset))
  {
    temp <- 0
    if (modelset[i] != "yule2rate" && modelset[i] != "rbvd")
        method <- parse(text = paste("I", modelset[i], "(x)", sep = ""))
    else
        method <- parse(text = paste("I", modelset[i], "(x, ints)", sep = ""))

    temp <- suppressWarnings(eval(method))
    #temp <- eval(method)

    #cat(method, temp$LH, "\n")
    res$model[i] <- modelset[i]
    if (modelset[i] == "pureBirth") {res$params[i] = "r1"
      res$np[i] <- 1
      res$mtype[i] <- "RC"
      res$mabb[i] <- "pb"}
    if (modelset[i] == "bd") {res$params[i] = "r1, a"
      res$np[i] <- 2
      res$mtype[i] <- "RC"
      res$mabb[i] <- "bd"}
    if (modelset[i] == "DDX") {res$params[i] = "r1, X"
      res$np[i] <- 2
      res$mtype[i] <- "RV"
      res$mabb[i] <- "DDX"}
    if (modelset[i] == "DDL") {res$params[i] = "r1, k"
      res$np[i] <- 2
      res$mtype[i] <- "RV"
      res$mabb[i] <- "DDL"}
    if (modelset[i] == "yule2rate") {res$params[i] = "r1, r2, ts"
      res$np[i] <- 3
      res$mtype[i] <- "RV"
      res$mabb[i] <- "y2r8"}
    if (modelset[i] == "rvbd") {res$params[i] = "r1, r2, a, ts"
      res$np[i] <- 4
      res$mtype[i] <- "RV"
      res$mabb[i] <- "rvbd"}
    if (modelset[i] == "yule3rate"){res$params[i] = "r1, r2, r3, st1, st2"
      res$np[i] <- 5
      res$mtype[i] <- "RV"
      res$mabb[i] <- "y3r8"}
    #fill res with estimated parameter values

    res$LH[i] <- temp$LH
    res$r1[i] <- temp$r1
    if (!is.null(temp$r2))
      res$r2[i] <- temp$r2
    else
      res$r2[i] <- NA
    if (!is.null(temp$a))
      res$a[i] <- temp$a
    else
      res$a[i] <- NA
    if (!is.null(temp$xp))
      res$xp[i] <- temp$xp
    else
      res$xp[i] <- NA
    if (!is.null(temp$k))
      res$k[i] <- temp$k
    else
      res$k[i] <- NA
    if (!is.null(temp$st))
      res$st[i] <- temp$st
    else
      res$st[i] <- NA
    if (!is.null(res$LH[i]))
      res$AIC[i] <- (-2*res$LH[i]) + 2*res$np[i]
    else
      res$AIC[i] <- NA
    if (!is.null(temp$st2))
      res$st2[i] <- temp$st2
    else
      res$st2[i] <- NA
    if (!is.null(temp$r3))
      res$r3[i] <- temp$r3
    else
      res$r3[i] <- NA
  }
  aicmin <- min(res$AIC)
  res$dAIC <- res$AIC - aicmin
  rcbest <- min(res$AIC[res$mtype == "RC"])
  rvbest <- min(res$AIC[res$mtype == "RV"])
  daicstat <- rcbest - rvbest

  summ <- list()
  summ$dAICrc <- daicstat
  summ$rcbest <- rcbest
  summ$rvbest <- rvbest
  summ$bestmodel <- res$mabb[res$AIC == min(res$AIC)]
  if (res$mabb[res$AIC == min(res$AIC)] == "DDX" && res$x[res$AIC == min(res$AIC)] < 0)
    summ$sp <- "RI"
  else if (res$mabb[res$AIC == min(res$AIC)] == "y2r8" && res$r1[res$AIC == min(res$AIC)] < res$r2[res$AIC == min(res$AIC)])
    summ$sp <- "RI"
  else if (res$mabb[res$AIC == min(res$AIC)] == "rvbd" && res$r1[res$AIC == min(res$AIC)] < res$r2[res$AIC == min(res$AIC)])
    summ$sp <- "RI"
  else if (res$mabb[res$AIC == min(res$AIC)] == "y3r8" && res$r1[res$AIC == min(res$AIC)] < res$r3[res$AIC == min(res$AIC)])
    summ$sp <- "RI"
  else if (res$mabb[res$AIC == min(res$AIC)] == "pb" || res$mabb[res$AIC == min(res$AIC)] == "bd")
    summ$sp <- "NC"
  else
    summ$sp <- "RD"

  for (i in 1:length(modelset))
  {
    eval(parse(text = paste("summ$L", res$model[i], " <- ", res$LH[i], sep = "")))
  }
  summ <- as.data.frame(summ)

  pvec <- c("LH", "AIC", "r1", "r2", "r3", "a", "x", "k", "st", "st2")

  if (verbose2 == TRUE)
  {
    sink(file, append = TRUE)
    cat(paste(summ$bestmodel), "\t", summ$dAICrc, "\t", paste(summ$sp), "\t")
    for (i in 1:length(modelset))
    {
      for (j in 1:length(pvec))
      {
        cat(eval(parse(text = paste("res$", pvec[j], "[i]", sep = ""))), "\t")
      }
    }
    cat("\n")
    sink()
  }

  return(summ)

}
