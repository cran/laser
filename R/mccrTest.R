`mccrTest` <-
function(CladeSize, NumberMissing, NumberOfReps, ObservedGamma=NULL, fname=NULL){
	
	if (is.null(fname)){
		v <- vector('list', NumberOfReps);
		for (i in 1:NumberOfReps)
			v[[i]] <- birthdeath.tree(b=.1,d=0, taxa.stop=CladeSize);	
		v <- lapply(v, old2new.phylo);	
		v <- lapply(v, drop.tip, as.character(sample(1:CladeSize, NumberMissing)));	
	}else{
		v <- read.tree(file = fname);	
	}
	NullGamma <- unlist(lapply(v, gammaStat));
	CriticalValue <- as.numeric(quantile(NullGamma, .05));
	if (!is.null(ObservedGamma)){
		p.value <- (length(NullGamma[NullGamma <= ObservedGamma]) + 1)/(1+ length(NullGamma));
		return(list(null.gamma = NullGamma, critical.value = CriticalValue, p.value=p.value));
	}else{
		return(list(null.gamma = NullGamma, critical.value = CriticalValue));
	}
	
}

