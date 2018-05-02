#===============================================================================================#
#' Random Walk Metropolis-Hastings Algorithm
#'
#' \code{RWMH} computes random draws of parameters using a specified proposal distribution.
#' The default is the normal distribution
#'
#' @param data data required for the posterior distribution
#' @param propob a list of mean and variance-covariance of the normal proposal distribution
#' (default: NULL i.e. internally generated)
#' @param posterior the posterior distribution. It is set to null in order to use the logit posterior.
#' The user can specify log posterior as a function of parameters and data (pars,data)
#' @param iter number of random draws desired (default: 15000)
#' @param burn burn-in period for the Random Walk MH algorithm (default: 1000)
#' @param start starting values of parameters for the MH algorithm.
#' It is automatically generated but the user can also specify.
#' @param prior the prior distribution (default: "Normal", alternative: "Uniform")
#' @param mu the mean of the normal prior distribution (default:0)
#' @param sig the variance of the normal prior distribution (default:10)
#' @return val a list of matrix of draws pardraws and the acceptance rate
#'
#' @examples
#' # data = genmle::dat_mroz
#' # propob<- lapl_aprx(data[,1],data[,-1])
#' # RWMHob<- RWMH(data=data,propob) # prior="Normal"
#' # system.time(RWMHob<- RWMH(data=data,propob)) # time operation
#' # RWMHob<- RWMH(data=data,prior="Uniform")
#' # RWMHob<- RWMH(data=data,iter=5000) # a higher number of draws
#'
#' @export

RWMH<- function(data,propob=NULL,posterior=NULL,iter=15000,burn=1000,start=NULL,prior="Normal",mu=0,sig=10){
  if(is.null(posterior)){
    logpost<- function(start,data) posterior(start,data,Log=T,mu=mu,sig=sig,prior=prior)
  #define posterior distribution
  }
  if(is.null(propob)){
    parrs = lapl_aprx(data[,1],data[,-1])
    propob=lapl_aprx2(parrs$mode,logpost,data=data) #approximate the actual posterior distribution
  }
  varprop = 1.5*propob$var
  npar = length(propob$mode)
  Mat = array(0, c(iter, npar))
  if(is.null(start)){
    start = MASS::mvrnorm(n=1,propob$mode,varprop)
  }
  e = 0.000001 # specify step unif(-e,e) e = 0.000001
  Mat[1,] = start; AccptRate<-0
  for(i in 2:iter){
      start= Mat[i-1,]
      prop = MASS::mvrnorm(n=1,start,varprop) + stats::runif(1,-e,e)#make a draw from proposal dist
      lpa = logpost(prop,data); lpb = logpost(start,data)
      accprob = exp(lpa-lpb)
      # the other part cancels out because the normal distribution is symmetric
      if(stats::runif(1)< accprob){
        Mat[i,]=prop
        AccptRate<- AccptRate +1
      }else{
        Mat[i,]=start
      }
  }
  cat("RWMH algorithm successful\n")
  val = list(Matpram=Mat[-c(1:burn),],AcceptanceRate = AccptRate/iter)
  return(val)
  }



#===============================================================================================#
#' Independence Metropolis-Hastings Algorithm
#'
#' \code{IndepMH} computes random draws of parameters using a specified proposal distribution.
#' The default is the normal distribution
#'
#' @param data data required for the posterior distribution
#' @param propob a list of mean and variance-covariance of the normal proposal distribution (default:NULL)
#' @param posterior the posterior distribution. It is set to null in order to use the logit posterior.
#' The user can specify log posterior as a function of parameters and data (pars,data)
#' @param iter number of random draws desired (default: 15000)
#' @param burn burn-in period for the Random Walk MH algorithm (default: 1000)
#' @param start starting values of parameters for the MH algorithm.
#' It is automatically generated but the user can also specify.
#' @param prior the prior distribution (default: "Normal", alternative: "Uniform")
#' @param mu the mean of the normal prior distribution (default:0)
#' @param sig the variance of the normal prior distribution (default:10)
#' @return val a list of matrix of draws pardraws and the acceptance rate
#'
#' @examples
#' # data = genmle::dat_mroz
#' # propob<- lapl_aprx(data[,1],data[,-1])
#' # IndepMHob<- IndepMH(data=data,propob,iter=3000) # prior="Normal"
#' # system.time(IndepMHob<- IndepMH(data=data),iter=3000) # time operation
#' # IndepMHob<- IndepMH(data=data,propob,prior="Uniform",iter=3000)
#' # IndepMHob<- IndepMH(data=data,propob,iter=3000) # a higher number of draws
#'
#' @export

IndepMH<- function(data,propob=NULL,posterior=NULL,iter=15000,burn=1000,start=NULL,prior="Uniform",mu=0,sig=10){
  if(is.null(posterior)){
    logpost<- function(start,data) posterior(start,data,Log=T,mu=mu,sig=sig,prior=prior)
    #define posterior distribution
  }
  if(is.null(propob)){
    parrs = lapl_aprx(data[,1],data[,-1])
    propob=lapl_aprx2(parrs$mode,logpost,data=data) #approximate the actual posterior distribution
    zj = 1.5*parrs$var
    if(any(eigen(zj)$values<=0)){
      varprop=1.5*parrs$var
    }else{
      varprop=1.5*propob$var
    }
    
  }else{
    varprop = 1.5*propob$var 
    }
  
  npar = length(propob$mode)
  Mat = array(0, c(iter, npar))
  if(is.null(start)){
    start = MASS::mvrnorm(n=1,propob$mode,varprop)
  }
  Mat[1,] = start; AccptRate<-0
  for(i in 2:iter){
    start= Mat[i-1,]
    prop = MASS::mvrnorm(n=1,propob$mode,varprop)#make a draw from proposal dist
    lpa = logpost(prop,data); lpb = logpost(start,data)
    accprob = exp(lpa-lpb)
    # the other part cancels out because the normal distribution is symmetric
    if(stats::runif(1)< accprob){
      Mat[i,]=prop
      AccptRate<- AccptRate +1
    }else{
      Mat[i,]=start
    }
  }
  cat("IndepMH algorithm successful\n")
  val = list(Matpram=Mat[-c(1:burn),],AcceptanceRate = AccptRate/iter)
  return(val)
}
