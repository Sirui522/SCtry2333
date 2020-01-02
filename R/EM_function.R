#' @title EM Clustering
#' @description Implement EM algorithm to estimate the parameters and mixture proportions.
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @param x sample from GMM
#' @param N number of clusters
#' @param it_num number of iterations
#' @return the clustering parameters and proportions
#' @examples
#' \dontrun{
#' data_x = c(rnorm(500,0),rnorm(500,2));
#' M = EM.R(data_x, N = 2);M
#' }
#' @export

EM.R = function(x, N, it_num){
  n = length(x)
  m = N
  miu <- runif(m)
  alpha <- c(0.2,0.8)
  prob <- matrix(rep(0,n*m),ncol=m)
  for (step in 1:it_num){
    # E step
    for (j in 1:m){
      prob[,j]<- sapply(x,dnorm,miu[j],1)
    }
    sumprob <- rowSums(prob)
    prob<- prob/sumprob

    oldmiu <- miu
    oldalpha <- alpha

    # M step
    for (j in 1:m){
      p1 <- sum(prob[ ,j])
      p2 <- sum(prob[ ,j]*x)
      miu[j] <- p2/p1
      alpha[j] <- p1/n

    }
    print(paste("[",step,"/",it_num,"]"))

    epsilo <- 1e-4
    if (sum(abs(miu-oldmiu))<epsilo &
        sum(abs(alpha-oldalpha))<epsilo) break
  }

  B = cbind(alpha,miu)
  colnames(B) = c("pi","mu")
  return(B)

}
