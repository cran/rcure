#' Prints rcure object
#' @description Prints rcure object
#' @param x an object of rcure
#' @param Var if it is TRUE, the funciton returns standard error estimated by bootstrap method. If it is False, the function only returns estimators of coefficients. By default, Var = TRUE
#' @param ROC if TRUE, the function plots ROC. By default, ROC = TRUE
#' @param ... further arguments to be passed to the printrcure function
#' @references Cai, C., Zou, Y., Peng, Y., & Zhang, J. (2012). smcure: An R-Package for estimating semiparametric mixture cure models. Computer methods and programs in biomedicine, 108(3), 1255-1260
#' @importFrom plyr ddply .
#' @export

printrcure <-
function(x,Var=TRUE, ROC=TRUE, ...)
{
  #library(plyr)
  if (is.null(Var))
    Var = TRUE
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  cat("\nCure probability model:\n")
  if (Var) {
    b <- array(x$b, c(length(x$b), 4))
    rownames(b) <- x$bnm
    colnames(b) <- c("Estimate", "Std.Error", "Z value",
                     "Pr(>|Z|)")
    b[, 2] <- x$b_sd
    b[, 3] <- x$b_zvalue
    b[, 4] <- x$b_pvalue
  }
  if (!Var) {
    b <- array(x$b, c(length(x$b), 1))
    rownames(b) <- x$bnm
    colnames(b) <- "Estimate"
  }
  print(b)
  cat("\n")
  cat("\nFailure time distribution model:\n")
  if (Var) {
    beta <- array(x$beta, c(length(x$beta), 4))
    rownames(beta) <- x$betanm
    colnames(beta) <- c("Estimate", "Std.Error", "Z value",
                        "Pr(>|Z|)")
    beta[, 2] <- x$beta_sd
    beta[, 3] <- x$beta_zvalue
    beta[, 4] <- x$beta_pvalue
  }
  if (!Var) {
    beta <- array(x$beta, c(length(x$beta), 1))
    rownames(beta) <- x$betanm
    colnames(beta) <- "Estimate"
  }
  print(beta)
  cat("\n")
  cat("\n Evaluation Metrics\n")
  if (!Var){
    metric <- x$metric
    if(ROC){
      plot(1-x$sensep[,2], x$sensep[,1], type = "s", xlab = "1 - Specificity", ylab = "Sensitivity")
    }
  }
  if (Var){
    metric.boot <- sapply(x$eva.boot, function(x) x$metric)
    ci <- apply(metric.boot, 1, quantile, probs = c(0.025, 0.975))
    metric <- cbind(x$metric, t(ci))
    colnames(metric)[1] <- c("Estimate")
    if(ROC){
      sensep.boot <- lapply(x$eva.boot, function(x) x$sensep)
      sensep.res <- lapply(sensep.boot,
                           function( sensep ){
                             cutpoints <- seq(0,1-1e-5,length = 100)
                             cutSenspe(sensep, cutpoints)
                           }
      )
      type <- rep(1:100, length(sensep.res))
      sensep.res <- do.call(rbind, sensep.res)
      rownames(sensep.res) <- NULL
      sensep.res <- data.frame(type, sensep.res)
      sen.ci <- ddply(sensep.res, .(type), function(res) c( quantile(res$sen, c(0.025, 0.5, 0.975)), quantile(res$spe, 0.5 ) ) )

      plot(1 - sen.ci[,5], sen.ci[,3], type = "s", xlab = "1 - Specificity", ylab = "Sensitivity")
      lines(1 - sen.ci[,5], sen.ci[,2], lty = 2, type = "s")
      lines(1 - sen.ci[,5], sen.ci[,4], lty = 2, type = "s")
    }
  }
  print(metric)
  invisible(x)
}
