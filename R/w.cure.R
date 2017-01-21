w.cure <- function(pi, delta, surv){
  delta + (1-delta) * logit.inv( log(surv) + logit(pi) )
#   delta + (1-delta) * pi * surv / (1 - pi + pi * surv)
}









