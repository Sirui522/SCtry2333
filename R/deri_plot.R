#' @title Plot of directional derivative
#' @description Plot the directional derivative from the clustering estimator to dirac measures.
#' @param data_x sample from GMM
#' @param G the estimate of clustering
#' @param inf inf bound of the range for searching for parameters
#' @param sup sup bound of the range for searching for parameters
#' @return maximum value of gradient function in some range.
#' @examples
#' \dontrun{
#' data_x = c(rnorm(500,0),rnorm(500,2));
#' M = CN.R(data_x, it_num = 20, inf = -1.5, sup = 3.5);
#' M_deri = deri.plot(M);M_deri
#' }
#' @export
deri.plot = function(data_x, G, inf, sup){

  f1 = function(x,c){  ##Gaussian likelihood function f(x,c)
    out = exp(-0.5*(x-c)^2)/sqrt(2*pi)
    return(out)
  }
  f2 = function(x,G){  ##Gaussian likelihood function f(x,G), G is a discrete distribution with G[,1]=weight and G[,2]=mass points
    out = 0
    for(j in 1:length(G[,1])){
      out = out+G[j,1]*f1(x,G[j,2])
    }
    return(out)
  }
  d = function(c,G){  ##Formula(14), function of directional derivative from G to \delta_c
    out = 0
    for(j in 1:length(data_x)){
      out = out + f1(data_x[j],c)/f2(data_x[j],G)
    }
    out = out - length(data_x)
    return(out)
  }

  #####################

  inta = seq(inf, sup ,0.05)  ##grids in the interval to **search for local maxima** in Step1, Algorithm2

  ######################

  plot_d = rep(0,length(inta))
  for(i in 1:length(inta)){
    plot_d[i] = d(inta[i],G)
  }
  plot(inta,plot_d,type = "l",xlim = c(inf,sup),ylab = "d(θ,G)")
  box()
  axis(side = 1, at = c(inf:sup))
  abline(h = 0, lty = 2)
  return(max(plot_d))
}


