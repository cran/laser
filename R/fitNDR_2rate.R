`fitNDR_2rate` <-
function (phy, eps = 0, combined = TRUE, rate.decrease = FALSE, rbounds=c(.001, .5)) 
{
    phy$node.label <- NULL
    root <- max(phy$edge) - phy$Nnode + 1
    node.list <- 1:max(phy$edge)
    node.list <- node.list[node.list != root]
    x <- branching.times(phy)
    res <- list()
    for (i in 1:length(node.list)) {
        r1 <- NULL
        r2 <- NULL
        z1 <- NULL
        z2 <- NULL
        z <- NULL
        z <- splitEdgeMatrix(phy, node.list[i])
        z1 <- z[z[, 6] == 1, ]
        z2 <- z[z[, 6] == 2, ]
        r1 <- getLambda.internal(z1, root, eps = eps, combined = combined, rbounds=rbounds)
        if (rate.decrease == FALSE) {
            r2 <- getLambda.internal(z2, root, rbounds=rbounds, eps = eps, 
                combined = combined)
        }
        else {
            r2 <- getLambda.internal(z2, root, rbounds = c(rbounds[1], r1$r), eps = eps, 
                combined = combined);
        }
        res$node[i] <- node.list[i]
        res$LH[i] <- r1$LH + r2$LH
        res$aic[i] <- (-2 * res$LH[i]) + 6
        res$r.1[i] <- r1$r
        res$lambda.1[i] <- r1$lambda
        res$LH.1[i] <- r1$LH
        res$r.2[i] <- r2$r
        res$lambda.2[i] <- r2$lambda
        res$LH.2[i] <- r2$LH
        res$eps[i] <- eps
    }
    res <- as.data.frame(res)
    return(res)
}

