`Iyule2rate` <-
function(x, ints = NULL)
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

