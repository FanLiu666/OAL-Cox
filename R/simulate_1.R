## Writen for R version 4.3.3

############################ Simulate data function ############################
# define the function note
#' simulate data
#' 
#' @param mean_x: mean of the covariates
#' @param sig_x: standard deviation of covariates
#' @param rho: correlation of covariates
#' @param n: sample size
#' @param p: number of covariates
#' @param lambda: weibull model parameters
#' @param alpha: weibull model parameters
#' @param theta: censored distribution parameters
#' @return list of datasets and variables

simulate_data <- function(mean_x, sig_x, rho, n, p, lambda, alpha, theta){
  # set covariates
  pC = pI = pP <- 2
  pS <- p - (pC + pI + pP)
  var.list <- c(paste("Xc", 1:pC, sep = ""), 
                paste("Xp", 1:pP, sep = ""), 
                paste("Xi", 1:pI, sep = ""), 
                paste("Xs", 1:pS, sep = ""))
  # set strength of relationship between covariates and outcome
  beta_v <-  c( 0.6, 0.6, 0.6, 0.6, 0, 0, rep(0, p-6) )
  # set strength of relationship between covariates and treatment
  alpha_v <- c( 1.0, 1.0, 0, 0, 1, 1, rep(0, p-6) )
  names(beta_v) = names(alpha_v) <- var.list
  # set coefficient of A
  bA <- 0.6
  
  ### define function for generating data
  expit <- function(x){ 
    pr <- ( exp(x) / (1 + exp(x)) ) 
    return(pr)
  }
  
  ### simulate data
  # simulate covariates
  Sigma_x <- matrix(rho*sig_x^2, nrow = length(var.list), ncol = length(var.list)) 
  diag(Sigma_x) <- sig_x^2
  Mean_x <- rep(mean_x, length(var.list))
  Data <- as.data.frame(mvrnorm(n = n, mu = Mean_x, Sigma = Sigma_x, empirical = FALSE)) 
  names(Data) <- var.list
  
  # simulate treatment A 
  gA_x <- rowSums(Data[ , var.list]*matrix(alpha_v, nrow = n, ncol = length(var.list), byrow = TRUE))
  pA <- expit(gA_x)
  Data$A <- as.numeric(runif(n = length(pA)) < pA)
  
  # simulate survival time
  # set.seed(20240510)
  u <- runif(n, 0, 1)
  epsilon <- rnorm(n = n, sd = sig_x)
  gT_X¦Â <- rowSums(Data[ , var.list]*matrix(beta_v, nrow = n, ncol = length(var.list), byrow = TRUE)) 
  gT_X¦Â <- gT_X¦Â + epsilon
  gT_XA <- gT_X¦Â - Data$A*bA
  Data$Sur <- (diag(-log(u)/lambda) %*% exp(-gT_XA))^(1/alpha)
  
  # simulate cencor time
  CT <- runif(n, 0, theta)
  
  # get the observed time
  Data$time <- pmin(Data$Sur, CT)
  
  # simulate status
  Data$status <- as.numeric(Data$Sur <= CT) 
  Data <- Data %>%
    select(-Sur)
  # cencor_freq <- sum(Data$status == 0)/n
  
  # normalize covariates to have mean 0 and standard deviation 1
  temp.mean <- colMeans(Data[ , var.list])
  Temp.mean <- matrix(temp.mean, ncol = length(var.list), nrow = nrow(Data), byrow = TRUE)
  Data[,var.list] <- Data[ , var.list] - Temp.mean
  temp.sd <- apply(Data[var.list], FUN = sd, MARGIN = 2)
  Temp.sd <- matrix(temp.sd, ncol = length(var.list), nrow = nrow(Data), byrow = TRUE)
  Data[var.list] <- Data[ , var.list] / Temp.sd
  rm(list = c("temp.mean", "Temp.mean", "temp.sd", "Temp.sd"))
  
  return(list(Data, var.list))
}