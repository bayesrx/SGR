##################################
###     SimBaS   Function      ###
##################################

#' SimBaS
#'
#' Implementation of SimBaS steps to find out which
#' part of the functions are bounded away from zero.
#'
#' @param f_mcmc_sample MCMC samples of the function f (n X n_mcmc)
#' @param Spatial_locations Spatial coordinates nX2
#'
#' @return a list containing:
#'         P_SimBaS <- local edge probability for each spatial locations.
#'         f_non_zero <- functional values where it is bounded away from zero.
#'         significant_locations <- locations where the function is away from zero.
#'
#'
SimBaS <- function(f_mcmc_samples,spatial_locations){

### Description: Implementation of SimBaS steps to find out which
###              part of the functions are bounded away from zero.

## Inputs: f_mcmc_sample <- MCMC samples of the function f (n X n_mcmc)
##         Spatial_locations: Spatial coordinates nX2

## Output: P_SimBaS <- local edge probability for each spatial locations.
##         f_non_zero <- functional values where it is bounded away from zero.
##         significant_locations <- locations where the function is away from zero.

n_mcmc <- dim(f_mcmc_samples)[2]

fS_MCMC_samples <- f_mcmc_samples
fS_mean <- apply(fS_MCMC_samples, 1, mean)
fS_std <- apply(fS_MCMC_samples, 1, sd)
Z <- abs(fS_MCMC_samples - matrix(rep(fS_mean,n_mcmc),nrow = length(fS_mean),ncol = n_mcmc))
Z <- diag(1/fS_std)%*%Z
Z <- apply(Z,2,max)

Z2 <- matrix(rep(Z,length(fS_mean)),nrow = length(fS_mean),ncol = n_mcmc,byrow = TRUE)
Z1 <- matrix(rep(abs(fS_mean/fS_std),n_mcmc),nrow = length(fS_mean),ncol = n_mcmc)
P_SimBaS <- ifelse(Z1 < Z2, 1, 0 )
P_SimBaS <- apply(P_SimBaS, 1, mean)
S_pos <- which(P_SimBaS < 0.05)

SimBaS_return <- list(P_SimBaS = P_SimBaS, f_non_zero = fS_mean[S_pos],
                      significant_locations = spatial_locations[S_pos,]) #S_pos=S_pos
return(SimBaS_return)

}

