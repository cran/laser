`getTipdata` <-
function (tipdata, phy)
{
    if (is.data.frame(tipdata)){
    	x <- as.vector(tipdata[,1]);
    	names(x) <- row.names(tipdata);	
    }else{
    	x<- tipdata;	
    }
    if (class(phy) != "phylo")
        stop("object \"phy\" is not of class \"phylo\"")
    if (is.null(phy$edge.length))
        stop("your tree has no branch lengths: invalid input")
    tmp <- phy$edge;
    nb.tip <- length(phy$tip.label)
    if (phy$Nnode != nb.tip - 1)
        stop("\"phy\" is not fully dichotomous")
    if (length(x) != nb.tip)
        stop("length of phenotypic and of phylogenetic data do not match")
    if (any(is.na(x)))
        stop("method can't be used with missing data...")
    phenotype <- as.numeric(rep(NA, nb.tip + phy$Nnode))
    names(phenotype) <- 1:max(phy$edge);
    if (is.null(names(x))) {
        phenotype[1:nb.tip] <- x;
    }
    else {
        if (!any(is.na(match(names(x), phy$tip.label))))
            phenotype[1:nb.tip] <- x[phy$tip.label]
        else {
            phenotype[1:nb.tip] <- x
            warning("the names of argument \"x\" and the names of the tip labels\ndid not match\n")
        }
    }
    phy$phenotype<-rep(NA, (nb.tip+phy$Nnode - 1))
    for (i in 1:length(phenotype))
    {
      phy$phenotype[as.numeric(phy$edge[,2])==names(phenotype[i])] <- phenotype[i]
    }
    return(phy);
}

