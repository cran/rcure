cutSenspe <- function(res, cut.point){
  cut.point <- quantile(res[,"cut"], cut.point)
  res <- res[order(res[,"cut"]),]
  res <- sapply(cut.point, function(cut.point){
    ind <- which(res[,"cut"] > cut.point)[1]
    colMeans( res[ind+ c(-1,0), ] )
  })
  t(res)
}
