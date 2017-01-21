auc.cure <- function(odds, risk, surv, delta){
  odds <- as.numeric(odds)
  risk <- as.numeric(risk)
  surv <- as.numeric(surv)
  delta <- as.numeric(delta)

  pi <- logit.inv(odds)
  w.est <- w.cure(pi, delta, surv)
  r1 <- outer(odds, odds, ">") + 0.5 * outer(odds, odds, "==")
  w1 <- outer(w.est, w.est, function(x,y) x * (1-y))
  num <- r1 * w1
  diag(num) <- 0
  den <- w1
  diag(den) <- 0

  sum(num) / sum(den)
}
