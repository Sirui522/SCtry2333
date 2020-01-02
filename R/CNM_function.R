#' @title CN Clustering
#' @description Implement constrained Newton method to estimate the parameters and mixture proportions.
#' @param data_x sample from GMM
#' @param it_num number of iterations
#' @param inf inf bound of the range for searching for parameters
#' @param sup sup bound of the range for searching for parameters
#' @return the clustering parameters and proportions
#' @examples
#' \dontrun{
#' data_x = c(rnorm(500,0),rnorm(500,2));
#' M = CN.R(data_x, it_num = 20, inf = -1.5, sup = 3.5);M
#' }
#' @export
CN.R = function(data_x, it_num, inf, sup){

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
  dd = function(c,G){  ##derivative of function:[ d ](above)
    out = 0
    for(j in 1:length(data_x)){
      out = out + (1/sqrt(2*pi))*(data_x[j]-c)*exp(-0.5*(data_x[j]-c)^2)/f2(data_x[j],G)
      }
    return(out)
  }
  delete = function(index,a){  ##function of deleting the same element as "a" in the vector "index"
    d = which(index == a)
    out = index[-d]
    return(out)
  }

  #####################

  G0 = matrix(c(rep(0.5,2),c(inf,sup)),2,2)  ##initial guess of the distribution with G0[,1]=weight and G0[,2]=mass points
  r = length(data_x)/1e6  ##parameter \gammma in formula (12)
  inta = seq(inf, sup ,0.05)  ##grids in the interval to **search for local maxima** in Step1, Algorithm2
  tol = 1e-5
  tol2 = 1e-3
  tol3 = 5e-3

  ######################

  G = G0  ##the estimate of mixing distribution with G[,1]=weight and G[,2]=mass points
  for(count in 1:it_num){

    ##Step 1
    ###Here, "theta" means mass points of the distribution and "pi" means the weight

    theta_plus = G[,2]  ##initialize support points of G_count

    record_d = rep(0,length(inta))
    for(i in 1:length(inta)){
      record_d[i] = d(inta[i],G)
    }  ##compute the function:[  d  ] for searching local maxima

    deri_d = rep(0,length(inta))
    for(i in 1:length(inta)){
      deri_d[i] = dd(inta[i],G)
    }  ##compute the derivative of function:[  d  ] for searching local maxima


    for(i in 1:(length(inta)-1)){ ##search for local maxima
      if(deri_d[i]>0 & deri_d[i+1]<0) {
        if((record_d[i]+record_d[i+1])/2>tol2){  ##if the d(theta_i,G) is very small, then no need to add it
          theta_plus = c(theta_plus,(inta[i]+inta[i+1])/2)  ##add the central point in each interval instead of the real local maxima
        }
      }
    }
    if(abs(record_d[length(inta)]-max(record_d))<tol) theta_plus = c(theta_plus,inta[length(inta)])
    if(abs(record_d[1]-max(record_d))<tol) theta_plus = c(theta_plus,inta[1])  ##the local maxima can appear at the boundary of the whole interval we search

    #if(count%%5==0) print(paste("theta_plus is:",theta_plus))



    ##Step 2
    pi_plus = c(G[,1],rep(0,length(theta_plus)-length(G[,2])))
    G_Plus1 = cbind(pi_plus,theta_plus)

    ###One step of the CN method, which is also the Step 2 of Algorithm 1 [Wang,2007], linked with Formula (8)-(10) and (12)
    S = matrix(data = 0, nrow = length(data_x)+1, ncol = length(theta_plus))
    deno = rep(0,length(data_x))
    for(i in 1:length(data_x)){
      deno[i] = f2(data_x[i],G_Plus1)
      for(j in 1:length(theta_plus)){
        S[i,j] = f1(data_x[i],theta_plus[j])/deno[i]
      }
    }  ##Formula (5)

    ###Formula (12) and the NNLS algorithm (Appendix A, [Wang, 2007])
    library(nnls)
    S[length(data_x)+1,] = rep(r^(-0.5),length(theta_plus))
    yb = c(rep(2,length(data_x)),r^(-0.5))
    pi_minus = nnls(S,yb)$x
    pi_minus = pi_minus/sum(pi_minus)  ##normalization



    ##Step 3
    index = c(1:length(pi_minus))
    for(j in c(1:length(pi_minus))){
      if(pi_minus[j]<tol3) index = delete(index,j)
    }

    theta_plus_new = theta_plus[index]
    pi_plus_new = pi_minus[index]
    G = cbind(pi_plus_new,theta_plus_new)

    print(paste("[",count,"/",it_num,"]"))

  }
  colnames(G) = c("pi","Mu")
  return(G)
}

