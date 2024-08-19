#' @title Outcome-adaptive Lasso for Cox Proportional Hazard Model
#'
#' @param data: survival data
#' @param A: treatment
#' @param time: observed survival time
#' @param status: survival status
#' @param covariates: confounders
#' @param lambda_vec: regularization parameter for OAL-Cox model
#' @param gamma_convergence_factor: converge factor
#' @param tau: time point of interest
#'
#' @return RASCE of OAL-Cox model, the optimal lambda_n and propensity score coefficients corresponding to this value
#' @export
#'
#' @examples

OAL_Cox <- function(data, A, time, status, covariates, lambda_vec, gamma_convergence_factor, tau){
  # sample number
  n <- nrow(data)
  # set vector of possible lambda's to try
  names(lambda_vec) <- as.character(lambda_vec)
  # lambda_n (n)^(gamma/2 - 1) = n^(gamma_convergence_factor)
  # get the gamma value for each value in the lambda vector that corresponds to convergence factor
  gamma_vals <- 2*( gamma_convergence_factor - lambda_vec + 1 )
  names(gamma_vals) <- names(lambda_vec)

  ### define some functions for weight calculation, and the wAMD
  create_weights <- function(fp, fA, fw){
    fw <- (fp)^(-1)
    fw[fA == 0] <- (1 - fp[fA == 0])^(-1)
    return(fw)
  }

  wAMD_function <- function(DataM, varlist, trt.var, wgt, beta){
    trt = untrt = diff_vec <- rep(NA,length(beta))
    names(trt) = names(untrt) = names(diff_vec) <- varlist
    for(jj in 1:length(varlist)){
      this.var <- paste("w", varlist[jj], sep = "")
      DataM[ , this.var] <- DataM[ , varlist[jj]] * DataM[ , wgt]
      trt[jj] <- sum(DataM[DataM[ , trt.var] == 1, this.var]) / sum(DataM[DataM[ , trt.var] == 1, wgt])
      untrt[jj] <- sum(DataM[DataM[ , trt.var] == 0, this.var]) / sum(DataM[DataM[ , trt.var] == 0, wgt])
      diff_vec[jj] <- abs(trt[jj] - untrt[jj])
    }
    wdiff_vec <- diff_vec * abs(beta)
    wAMD <- c(sum(wdiff_vec))
    ret <- list(diff_vec = diff_vec, wdiff_vec = wdiff_vec, wAMD = wAMD)

    return(ret)
  }

  # estimate Cox PH model
  y.form <- formula(paste("Surv(time, status) ~ A + ", paste(covariates, collapse = " + ")))
  lm.Y <- coxph(y.form, data = data)
  betaXY <- coef(lm.Y)[covariates]

  ### Want to save wAMD and propensity score coefficients for each lambda value
  wAMD_vec <- rep(NA, length(lambda_vec))
  names(wAMD_vec) <- names(lambda_vec)
  coeff_XA <- as.data.frame(matrix(NA, nrow = length(covariates), ncol = length(lambda_vec)))
  names(coeff_XA) <- names(lambda_vec)
  rownames(coeff_XA) <- covariates

  ### Run OAL-Cox model for each lambda value
  # weight model with all possible covariates included, this is passed into lasso function
  w.full.form <- formula(paste("A ~ ", paste(covariates, collapse = " + ")))
  for( lil in names(lambda_vec) ){
    il <- lambda_vec[lil]
    ig <- gamma_vals[lil]

    # create the OAL-Cox penalty with coefficient specific weights determined by Cox PH model
    OAL_Cox_pen <- adaptive.lasso(lambda = n^(il), al.weights = abs(betaXY)^(-ig) )
    # run OAL-Cox model with appropriate penalty
    logit_OAL_Cox <- lqa.formula(w.full.form,
                                 data = data,
                                 penalty = OAL_Cox_pen,
                                 family = binomial(logit))

    # generate propensity score
    data[ , paste("f.pA", lil, sep = "")] <- predict.lqa(logit_OAL_Cox)$mu.new
    # save propensity score coefficients
    coeff_XA[covariates, lil] <- coef(logit_OAL_Cox)[covariates]

    # create inverse probability of treatment weights
    data[ , paste("w", lil, sep = "")] <- create_weights(fp = data[ , paste("f.pA", lil, sep = "")], fA = A)

    # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
    wAMD_vec[lil] <- wAMD_function(DataM = data,
                                   varlist = covariates,
                                   trt.var = "A",
                                   wgt = paste("w", lil, sep = ""),
                                   beta = betaXY)$wAMD
  } # close loop through lambda values

  # find the lambda value that creates the smallest wAMD
  tt = which.min(wAMD_vec)

  # print out the coefficients for the propensity score that corresponds with smallest wAMD value to select covariates
  coeff <- coeff_XA[tt]
  coeff_OAL_Cox <- ifelse(abs(coeff) <= 1e-8, 0, 1)
  select_var <- covariates[which(coeff_OAL_Cox != 0)]

  # print out RASCE corresponding to the covariates selected by OAL-Cox
  RASCE_OAL_Cox <- calc_ipw_rmst(data = data,
                                 A = A,
                                 time = time,
                                 status = status,
                                 covariates = select_var,
                                 tau = tau)

  return(list(RASCE_OAL_Cox, tt, coeff_OAL_Cox, data))
}
