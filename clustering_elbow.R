library('amap')

get_inter_cluster_dists <- function(centers) {
  toRet <- 0
  for(i in 1:nrow(centers)) {
    for(j in i:nrow(centers)) {
      toRet <- toRet + (dist(rbind(centers[i,], centers[j,]))^2)
    }
  }
  return(mean(toRet) / length(toRet))
}

elbow_scores <- c()
for(i in 1:10) {
  tmp <- Kmeans(t(dd), centers = i, method = "maximum", iter.max = 1000)
  inter_cluster_dists <- get_inter_cluster_dists(tmp$centers)
  elbow_scores <- c(elbow_scores, sum(tmp$withinss) / inter_cluster_dists )
}
plot(1:10, elbow_scores, type='l', col='blue',
     xlab = 'Number of centers', ylab = 'Score',
     main = 'Score with increased # centers in kmeans')
