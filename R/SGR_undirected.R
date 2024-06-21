## Spatial Graphical Regression function (undirected setting)

SGR_undirected <- function(features,spatial_locations,ncores,N_MCMC,N_BURNIN,N_thin){

### Description: Spatial Graphical Regression for undirected setting
###              This code is based on MCMC algorithm

## Inputs:
## Mandatory inputs
## (1) features: This is gene expression matrix of dimension nXp
## (2) Spatial_locations: Spatial coordinates nX2
## (3) ncores: Number of cores
## (4) N_MCMC: Number of MCMC iterations to be saved
## (5) N_BURNIN: Number of burnin iterations (before start to save)
## (6) N_thin: Number of iterations to skip between two saved iterations

## Outputs:
## A list containing mcmc samples of for a node specific regression.
## (1) Beta_post: Posterior samples of beta
## (2) gamma_prob_post: Posterior samples of probability of selection indicators
## (3) sigma_epsilon_post: Posterior samples of error variance

#######################
###    libraries    ###
#######################
# library(fastBayesReg)
# library(BayesGPfit)
# library(ggplot2)
# library(reshape2)
# library(rlist)
# library(foreach)
# library(iterators)
# library(parallelly)
# library(doParallel)

################################################################################
  ## Spatial Graphical Regression function (single regression)

#' Spatial Graphical Regression function (undirected setting)
#'
#' Spatial Graphical Regression for undirected setting
#'      This code is based on MCMC algorithm
#'
#' @param features This is gene expression matrix of dimension nXp
#' @param Spatial_locations Spatial coordinates nX2
#' @param ncores Number of cores
#' @param N_MCMC Number of MCMC iterations to be saved
#' @param N_BURNIN Number of burnin iterations (before start to save)
#' @param N_thin Number of iterations to skip between two saved iterations
#'
#' @return
#' A list containing mcmc samples of for a node specific regression.
#' (1) Beta_post: Posterior samples of beta
#' (2) gamma_prob_post: Posterior samples of probability of selection indicators
#' (3) sigma_epsilon_post: Posterior samples of error variance
#'
#'
#'
#'

  SGR_undirected_singleReg <- function(features,spatial_locations,node,N_MCMC,N_BURNIN,N_thin){

    ### Description: Spatial Graphical Regression for single regression with undirected setting
    ###              This code is based on MCMC algorithm

    ## Inputs:
    ## Mandatory inputs
    ## (1) features: This is gene expression matrix of dimension nXp
    ## (2) Spatial_locations: Spatial coordinates nX2
    ## (3) node: Mention the node on the regression needs to be performed
    ## (4) N_MCMC: Number of MCMC iterations to be saved
    ## (5) N_BURNIN: Number of burnin iterations (before start to save)
    ## (6) N_thin: Number of iterations to skip between two saved iterations

    ## Outputs:
    ## A list containing mcmc samples of for a node specific regression.
    ## (1) Beta_post: Posterior samples of beta
    ## (2) gamma_prob_post: Posterior samples of probability of selection indicators
    ## (3) sigma_epsilon_post: Posterior samples of error variance

    #######################
    ###    libraries    ###
    #######################
    # library(fastBayesReg)
    # library(BayesGPfit)
    # library(ggplot2)
    # library(reshape2)
    # library(rlist)

    Y <- features  ## Count matrix (n x p)
    S <- spatial_locations  ##location matrix (n x 2)

    n <- dim(features)[1]  ## number of locations
    p <- dim(features)[2]  ## number of genes

    ################################################################
    ######          Basis expansion for Beta function          #####
    ################################################################
    ### Generate eigen functions
    a_beta <- 0.01
    b_beta <- 10
    poly_degree_S <- 25

    ## eigen function for spatial function
    Psi_beta_S <- GP.eigen.funcs.fast(S, poly_degree = poly_degree_S, a = a_beta, b = b_beta) #Eigen function matrix (n x L_S)
    L_S <- ncol(Psi_beta_S)

    ### eigen values for spatial function
    lambda_beta_S <- GP.eigen.value(poly_degree = poly_degree_S, a = a_beta, b = b_beta, d = ncol(S)) #Eigen values vector of length L_S

    ### Bases for spatial eigen function
    Bases_beta_S <- t(Psi_beta_S)*sqrt(lambda_beta_S) #Basis matrix L_S x n
    ########################        XXXXXXXXXXXXXXXXX        ########################

    ################################################################
    ######          Basis expansion for eta function           #####
    ################################################################
    ### Generate eigen functions
    a_eta <- runif(1,0.1,0.25)
    b_eta <- runif(1,0.01,0.1)
    poly_degree_S <- 25

    ## eigen function for spatial function
    Psi_eta_S <- GP.eigen.funcs.fast(S, poly_degree = poly_degree_S, a = a_eta, b = b_eta) #Eigen function matrix (n x L_S)
    L_S <- ncol(Psi_eta_S)

    ### eigen values for spatial function
    lambda_eta_S <- GP.eigen.value(poly_degree = poly_degree_S, a = a_eta, b = b_eta, d = ncol(S)) #Eigen values vector of length L_S

    ### Bases for spatial eigen function
    Bases_eta_S <- t(Psi_eta_S)*sqrt(lambda_eta_S) #Basis matrix L_S x n
    ########################        XXXXXXXXXXXXXXXXX        ########################

    ### Initialization of the selection parameter gamma
    index <- c(1:p);index <- index[-node]
    X_new <- t(Bases_beta_S)*as.vector(Y[,index[1]])
    index_new <- index[-c(index[1])]
    for(j in index_new){
      X_new <- cbind(X_new, t(Bases_beta_S)*as.vector(Y[,j]) )
    }
    X_new <- cbind(X_new,t(Bases_eta_S))
    Y_node <- Y[,node]

    fit <- fast_mfvb_normal_lm(Y_node,X_new)

    Beta_mat <- matrix(0,nrow=n,ncol=p)

    gamma <- rep(0,p)
    gamma_prob <- rep(0.5,p);gamma_prob[node]=0
    for (j in index) {
      if(j<node){
        aa <- 1+(j-1)*L_S;bb <- j*L_S
        Beta_mat[,j] <- crossprod(Bases_beta_S,fit$post_mean$betacoef[aa:bb])
        gamma[j] <- ifelse( sqrt(mean(Beta_mat[,j]^2))>0.1 ,1,0)
      }
      if(j==node)
      {
        Beta_mat[,j] <- 0
        gamma[j] <- 0
      }
      if(j>node){
        k <- j-1
        aa <- 1+(k-1)*L_S;bb <- k*L_S
        Beta_mat[,j] <- crossprod(Bases_beta_S,fit$post_mean$betacoef[aa:bb])
        gamma[j] <- ifelse( sqrt(mean(Beta_mat[,j]^2))>0.1 ,1,0)
      }
    }
    Noise_var <- fit$post_mean$sigma2_eps
    Prior_var <- fit$post_mean$tau2*Noise_var

    ### MCMC settings
    N_total <- N_MCMC + N_BURNIN
    N_effsamp <- (N_total - N_BURNIN)/N_thin

    ## output ##
    Beta_out <- array(0,c(n,p,N_effsamp))  ##Posterior sample matrix for Beta_tilde = Beta*gamma
    gamma_prob_out <- matrix(0,nrow=p,ncol=N_effsamp)  ##Posterior probability matrix of gamma
    sigmasq_eps_out <- rep(0,N_effsamp)  ##Posterior sample matrix of error variance

    ### Start of MCMC iteration
    for (m in 1:N_total) {

      if( sum(gamma==1)==0 ){
        if(m < 100){
          gamma <- rep(1,p);gamma[node]=0;
          #print("None of the covariate nodes are connected to the response node.")
        }
      }

      if( sum(gamma==1) > 1 ){
        pos <- which(gamma == 1)
        pos1<-pos[1]
        X_new <- t(Bases_beta_S)*as.vector(Y[,pos1])
        for (j in 2:length(pos)) {
          X_new <- cbind(X_new, t(Bases_beta_S)*as.vector(Y[,pos[j]]) )
        }
        X_new <- cbind(X_new,t(Bases_eta_S))
        fit <- fast_normal_lm(Y_node,X_new,mcmc_sample = 1,burnin = 0)
        Noise_var <- fit$post_mean$sigma2_eps
        Prior_var <- fit$post_mean$tau2*Noise_var

        for (j in 1:length(pos)) {
          aa <- 1+(j-1)*L_S;bb <- j*L_S
          Beta_mat[,pos[j]] <- crossprod(Bases_beta_S,fit$post_mean$betacoef[aa:bb])
        }

        aa_eta <- 1+(length(pos))*L_S;bb_eta <- (length(pos)+1)*L_S
        spatial_rm <- crossprod(Bases_eta_S,fit$post_mean$betacoef[aa_eta:bb_eta])

        pos0 <- which(gamma == 0)
        if(sum(pos0) > 0){
          u <- numeric(0)
          for (l in 1:L_S) {
            u[l] <- rnorm(1,0,sqrt(Prior_var*lambda_beta_S[l]))
          }
          beta_prior <- crossprod(Bases_beta_S, as.matrix(u) )

          for (j in 1:length(pos0)) {
            Beta_mat[,pos0[j]] <- beta_prior
          }
        }
        Beta_mat[,node] <- rep(0,n)

        ### Estimation of probability of selection indicators
        for (j in index) {
          gamma[j] <- 1
          c <- sum((Y_node - rowSums(Y*(Beta_mat%*%diag(gamma)))-spatial_rm)^2)
          c <- log(gamma_prob[j]+0.001) - (c/(2*fit$post_mean$sigma2_eps))  ##Added 0.001 for numerical stability

          gamma[j] <- 0
          d <- sum((Y[,node] - rowSums(Y*(Beta_mat%*%diag(gamma)))-spatial_rm)^2)
          d <- log(1 - gamma_prob[j]+0.001) - (d/(2*fit$post_mean$sigma2_eps))  ##Added 0.001 for numerical stability

          if(c==0 && d==0 ){gamma_prob[j] <- 0}
          else{ gamma_prob[j] <- 1/(1+exp(d-c)) }
        }
        gamma <- ifelse(gamma_prob>0.5,1,0)

      }

      if( sum(gamma==1) == 1 ){
        pos <- which(gamma==1)
        X_new <- t(Bases_beta_S)*as.vector(gamma[pos]*Y[,pos])
        X_new <- cbind(X_new,t(Bases_eta_S))
        fit <- fast_normal_lm(Y_node,X_new,mcmc_sample = 1,burnin = 0)
        Noise_var <- fit$post_mean$sigma2_eps
        Prior_var <- fit$post_mean$tau2*Noise_var

        Beta_mat[,pos] <- crossprod(Bases_beta_S,fit$post_mean$betacoef[1:L_S])

        aa_eta <- (1+L_S);bb_eta <- 2*L_S
        spatial_rm <- crossprod(Bases_eta_S,fit$post_mean$betacoef[aa_eta:bb_eta])

        pos0 <- which(gamma == 0)
        if(sum(pos0) > 0){
          u <- numeric(0)
          for (l in 1:L_S) {
            u[l] <- rnorm(1,0,sqrt(Prior_var*lambda_beta_S[l]))
          }
          beta_prior <- crossprod(Bases_beta_S, as.matrix(u) )

          for (j in 1:length(pos0)) {
            Beta_mat[,pos0[j]] <- beta_prior
          }
        }

        Beta_mat[,node] <- rep(0,n)

        ### Estimation of probability of selection indicators
        gamma_copy <- gamma
        gamma_copy[pos] <- 1
        c <- sum((Y[,node] - Y[,pos]*Beta_mat[,pos] - spatial_rm)^2)
        c <- log(gamma_prob[pos]+0.001) - (c/(2*fit$post_mean$sigma2_eps)) ##Added 0.001 for numerical stability

        gamma_copy[pos] <- 0
        d <- sum((Y[,node] - spatial_rm)^2)
        d <- log(1 - gamma_prob[pos]+0.001) - (d/(2*fit$post_mean$sigma2_eps)) ##Added 0.001 for numerical stability

        if(c==0 && d==0 ){gamma_prob[pos] <- 0}
        else{ gamma_prob[pos] <- 1/(1+exp(d-c)) }

        gamma[pos] <- ifelse(gamma_prob[pos] > 0.5,1,0)

      }

      if(m%%100 == 0) print(m)

      # -- save output -- #
      if(m > N_BURNIN && m%%N_thin == 0){

        Beta_out[,,(m-N_BURNIN)/N_thin] = Beta_mat
        gamma_prob_out[,(m-N_BURNIN)/N_thin] = gamma_prob
        sigmasq_eps_out[(m-N_BURNIN)/N_thin] = Noise_var

      }
    }

    SGR_undirected_singleReg_post <- list(list(Beta_post = Beta_out,gamma_prob_post=gamma_prob_out,sigmasq_eps_post=sigmasq_eps_out))
    return(SGR_undirected_singleReg_post)

  }
################################################################################


Y <- features  ## Count matrix (n x p)
S <- spatial_locations  ##location matrix (n x 2)

n <- dim(features)[1]  ## number of locations
p <- dim(features)[2]  ## number of genes

# make new cluster
cl = parallelly::makeClusterPSOCK(ncores, autoStop = TRUE)

# Register parallel backend
doParallel::registerDoParallel(cl)

# parallel for loop
estimates = foreach::foreach(
  j = 1:p,
  .packages = c('fastBayesReg', 'BayesGPfit', 'ggplot2',
                'reshape2', 'sparsevb', 'rlist', 'Matrix', 'matrixcalc'),
  .export = NULL,
  .combine = c) %dopar% {
  #.combine = list) %dopar% {

    Node_reg <- SGR_undirected_singleReg(Y,S,node=j,N_MCMC = 100,N_BURNIN = 100,N_thin = 1)

    Node_reg
  }

# close connections
on.exit(stopCluster(cl))

return(estimates)

}

