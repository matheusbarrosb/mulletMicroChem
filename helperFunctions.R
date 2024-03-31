### Auxiliary functions ####

getInits <- function (mu, sigma, Nchains) {

  initialMix <- function(chainID = 1, # helper function
                       mu, sigma,
                       N_chains) {
  
  out = list(mu = c((mu[1]), mu[2]),
           sigma = c(sigma[1], sigma[2]),
           chainID = chainID)
  return(out)
}
  inits <- lapply(1:Nchains,
                   function(id) initialMix(chainID = id,
                                           mu = c(0.9,1.3),
                                           sigma = c(0.3,0.01)))
  return(inits)
}































