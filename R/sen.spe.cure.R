#' @importFrom stats stepfun
sen.spe.cure <- function(odds, risk, surv, delta, cut.point = NULL, verbose = T){
  odds <- as.numeric(odds)
  risk <- as.numeric(risk)
  surv <- as.numeric(surv)
  delta <- as.numeric(delta)
  n <- length(odds)

  pi <- logit.inv(odds)
  w.est <- w.cure(pi, delta, surv)
  den <- sum(w.est)


  odds.o  <- sort(odds)
  w.est.o <- w.est[order(odds)]

  sen <- c(1, 1 - cumsum(w.est.o) / sum(w.est), 0)
  spe <- c(0,  cumsum(1 - w.est.o) / sum(1-w.est), 1 )
  odds.o <- c( min(odds.o) - 1e-5, odds.o, max(odds.o) + 1e-5)

  if(max(table(odds.o)) > 1){  # Handle ties
    if(verbose){ cat("there are ties in estimated odds\n") }
    sen <- unlist( tapply(sen, odds.o, function(x) x[c(length(x))]) )
    spe <- unlist( tapply(spe, odds.o, function(x) x[c(length(x))]) )
    a0 <- sum(diff(spe) * sen[-1] )
    a1 <- sum(diff(spe) * sen[-length(sen)] )
    auc0 <- mean( c(a0,a1) )
    sen[length(sen) - 1] <- sen[length(sen) - 2]
    sen[1] <- sen[2]
    sen <- c(1, sen)
    spe <- c(0, spe)

  }else{
    auc0 <- sum(diff(spe) * sen[-n] )
  }
  auc1 <- auc.cure(odds, risk, surv, delta) # Check the results is equivalent to directly caculation

  den <- ( sum(w.est) * sum(1-w.est) )
  rate <- ( den - sum(w.est * (1-w.est)) ) / den
  auc2 <- auc0 * rate

  if(verbose){
    cat(rate, "AUC is", round(auc0,3), "for integration with correction", round(auc2,3), ". AUC is ", round(auc1,3), "for directly caculation\n" )
    roc <- stepfun(rev( 1 - spe), c(rev(sen),1) )
    plot(roc, xlim = c(0,1), xlab = "1 - Specificity", ylab = "Sensitivity", main = "ROC curve" )
  }

  options(warn=-1)
  res <- cbind(sen, spe, cut = c(min(odds.o) - 1, unique(odds.o) ) )
  options(warn=0)

  if(! is.null(cut.point) ){
    res <- cutSenspe(res, cut.point)
  }
  res
}

