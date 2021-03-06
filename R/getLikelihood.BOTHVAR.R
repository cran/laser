`getLikelihood.BOTHVAR` <-
function(lam0, k, mu0, z, SpeciationTimes, Tmax)
{

  x <- SpeciationTimes;
  L2 <- 0; L3 <- 0; 
  L1<- sum(log(2:(length(x) -1))) + sum(log(lambdaFx(lam0, k, x[3:length(x)]))); 

  L2 <- (2*log(pTt.BOTHVAR(lam0, k, mu0, z, x[1], Tmax)) + 2*rhoFxn.BOTHVAR(lam0, k, mu0, z, x[1], Tmax));
  for (i in 3:length(x)){
  	L3 <- (L3 + 2*log(pTt.BOTHVAR(lam0, k, mu0, z, x[i], Tmax)) + rhoFxn.BOTHVAR(lam0, k, mu0, z, x[i], Tmax));
  }
  #cat(L1, L2, L3, sep='\n')
  Likelihood <- L1 + L2 + L3;
  Likelihood; 
}

