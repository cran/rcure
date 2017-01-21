#' Evaluate Mixture cure model fitted by rcure
#'
#' @export
#' @param time survival time
#' @param delta event status (1 = event)
#' @param X design matrix of survival part
#' @param beta parameters of survival part
#' @param Z design matrix of cure part
#' @param b parameters of cure part
#' @param model the type of survival model ("PH", "PO","Normal")
#' @param baseline whether the surv is baseline survival probability
#' @param surv survival probability
eva_cure <- function(time,delta,X,beta,Z,b, surv, model, baseline = T){
  est.risk <- X %*% beta
  est.odds <- Z %*% b
  est.pi   <- logit.inv(est.odds)
  if(model == "PH" & baseline == T){ est.surv <- surv ^ exp( est.risk) }else{ est.surv <- surv}
  est.w <- w.cure(est.pi, delta, est.surv)

  k1  <- k.ind(risk = est.risk, pi = est.pi, model)
  k1w <- k.ind(risk = est.risk, pi = est.w, model)
  c1  <- cindex(risk = est.risk, pi = est.pi, time, delta )
  c1w <- cindex(risk = est.risk, pi = est.w,  time, delta )
  a1  <- auc.cure(est.odds, est.risk, est.surv, delta)
  res1 <- sen.spe.cure(est.odds, est.risk, est.surv, delta, verbose = F )

  res <- list( para = c( beta, b),
               metric = c(AUC = a1, K = k1, C = c1w),
               sensep = res1 )
  res
}


