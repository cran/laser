`splitEdgeMatrix` <-
function (phy, node) 
{
    x <- branching.times(phy)
    rootnode <- length(phy$tip.label) + 1
    phy$tag <- rep(1, nrow(phy$edge))
    if (node >= rootnode) {
        node.desc <- node
        pos <- 1
        phy$tag[phy$edge[, 1] == node.desc[1]] <- 2
        while (pos != (length(node.desc) + 1)) {
            temp <- node.sons(phy, node.desc[pos])
            temp <- temp[temp > rootnode]
            for (k in 1:length(temp)) {
                phy$tag[phy$edge[, 1] == temp[k]] <- 2
            }
            node.desc <- c(node.desc, temp)
            pos <- pos + 1
        }
    }
    else if (node > 0) 
        phy$tag[phy$edge[, 2] == node] <- 2
    z <- cbind(phy$edge, gsr(phy$edge[, 1], names(x), x), phy$edge.length, 
        phy$phenotype, phy$tag)
    z <- matrix(as.numeric(z), dim(z))
    z <- as.data.frame(z)
    return(z)
}

