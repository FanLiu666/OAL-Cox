#' AL_Cox function
#' see "Zhang H H, Lu W. Adaptive Lasso for Cox's proportional hazards model[J]. Biometrika, 2007, 94(3): 691-703."
#'
#' @param data: survival data
#' @param A: treatment
#' @param time: observed survival time
#' @param status: survival status
#' @param covariates: confounders
#' @param lambda_vec: regularization parameter for AL-Cox model
#' @param tau: time point of interest
#' @return RASCE of AL-Cox model and coefficient of covariates

AL_Cox <- function(data, A, time, status, covariates, lambda_vec, tau){
  # sample number
  n <- nrow(data)

  # calculate lambda
  lambda = n^(lambda_vec)

  # calculate penalty
  y.form <- formula(paste("Surv(time, status) ~ A + ", paste(covariates, collapse = " + ")))
  Cox.Y <- coxph(y.form, data = data)
  betaXY <- coef(Cox.Y)
  penalty_factor <- 1/abs(betaXY)

  # select x
  data_x <- data %>%
    select(-time, -status)
  data_x <- data_x[, c("A", covariates)]
  x <- as.matrix(data_x)

  # select y
  time <- as.double(time)
  status <- as.double(status)
  y <- Surv(time, status)

  # run AL-Cox model
  AL_Cox <- cv.glmnet(x, y, alpha = 1, gamma = 1,
                      lambda = lambda,
                      family = "cox",
                      type.measure = "deviance",
                      penalty.factor = penalty_factor,
                      nfolds=10)
  coeff <- coef(AL_Cox, s =  AL_Cox$lambda.min)
  coeff <- coeff[-1, ]
  coeff_AL_Cox <- ifelse(abs(coeff) <= 1e-8, 0, 1)
  select_var <- names(coeff_AL_Cox)[which(coeff_AL_Cox != 0)]

  # calculate RASCE for AL-Cox method
  RASCE_AL_Cox <- calc_ipw_rmst(data = data,
                                A = A,
                                time = time,
                                status = status,
                                covariates = select_var,
                                tau = tau)

  return(list(RASCE_AL_Cox, coeff_AL_Cox))
}
