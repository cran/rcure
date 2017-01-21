
#' Evaluate Mixture cure model fitted by rcure
#'
#' @export
#' @param fit an extended smcure object
#' @param model the type of survival model ("PH", "PO","Normal")
eva_rcure <- function(fit, model){

  X <- fit$X
  Z <- fit$Z
  time <- fit$Time
  delta <- fit$status
  beta <- fit$beta
  b <- fit$b
  surv <- fit$s
  eva_cure(time,delta,X,beta,Z,b,surv,model)

}
