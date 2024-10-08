#' @title IPW-weighting Model
#'
#' @param data: survival data
#' @param A: treatment
#' @param time: observed survival time
#' @param status: survival status
#' @param covariates: confounders
#' @param tau: time point of interest
#'
#' @return IPW-RMST
#' @export
#'
#' @examples

calc_ipw_rmst <- function(data, A, time, status, covariates, tau){
  # calculate propensity score and IPW
  w.full.form <- formula(paste("A ~ ", paste(covariates, collapse = " + ")))
  logit <- glm(w.full.form, data=data, family=binomial(link='logit'))
  pred <- predict(logit, type='response')
  data$weight <- A/pred + (1-A)/(1-pred)

  # calculate IPW-RMST
  ipw_rmst <- akm_rmst(time=time,
                       status=status,
                       group=as.factor(A),
                       weight=data$weight,
                       tau = tau)$Est.
  ipw_rmst <- as.numeric(ipw_rmst)

  return(ipw_rmst)
}
