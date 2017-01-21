
cindex <- function(risk, pi, time, delta, tau = NULL){
  risk <- as.numeric(risk)
  pi   <- as.numeric(pi)
  n <- length(risk)
  t1 <- outer(time, time, "<"  )
  t2 <- rep(1, n)
  if(! is.null(tau) )t2 <- time < tau  # Boundry of time
  r1 <- outer(risk, risk, function(x,y) (x>y) + 0.5 * (x==y))
  p1 <- outer(pi, pi ,"*")
  num <- delta * t1 * t2 * r1 * p1
  diag(num) <- 0
  dev <- delta* t1 * t2 * p1
  diag(dev) <- 0
  sum(num) / sum(dev)
}





