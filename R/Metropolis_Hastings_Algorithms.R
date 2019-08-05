#===============================================================================================#
#' Random Walk Metropolis-Hastings Algorithm
#'
#' \code{RWMH} computes random draws of parameters using a specified proposal distribution.
#' The default is the normal distribution
#'
#' @param data data required for the posterior distribution. First column is the outcome
#' @param propob a list of mean and variance-covariance of the normal proposal distribution
#' (default: NULL i.e. internally generated). Write list as propob=list(mode=mode,var=variance-covariance)
#' @param posterior the posterior distribution. It is set to null in order to use the logit posterior.
#' The user can specify log posterior as a function of parameters and data (pars,data). For a more
#' flexible and generic implementation, use \code{\link{rwMHgen}}
#' @param iter number of random draws desired
#' @param burn burn-in period for the Random Walk MH algorithm
#' @param vscale a positive value to scale up or down the variance-covariance matrix in
#' the proposal distribution
#' @param start starting values of parameters for the MH algorithm.
#' It is automatically generated from the proposal distribution but the user can also specify.
#' @param prior the prior distribution (default: "Normal", alternative: "Uniform")
#' @param mu the mean of the normal prior distribution (default:0)
#' @param sig the variance of the normal prior distribution (default:10)
#' @return val a list of matrix of draws Matpram and the acceptance rate
#'
#' @examples
#' y = indicat(faithful$waiting,70)
#' x = scale(cbind(faithful$eruptions,faithful$eruptions^2))
#' data = data.frame(y,x); propob<- lapl_aprx(y,x)
#' RWMHob_n<- RWMH(data=data,propob,iter = 102, burn = 2) # prior="Normal"
#' RWMHob_u<- RWMH(data=data,propob,prior="Uniform",iter = 102, burn = 2)
#' par(mfrow=c(3,1));invisible(apply(RWMHob_n$Matpram,2,function(x)plot(density(x))))
#' invisible(apply(RWMHob_u$Matpram,2,function(x)plot(density(x))));par(mfrow=c(1,1))
#'
#' @export

RWMH<- function(data,propob=NULL,posterior=NULL,iter=1500,burn=500,vscale=1.5,
                start=NULL,prior="Normal",mu=0,sig=10){
  if(is.null(posterior)){
    logpost<- function(start,data) posterior(start,data,Log=TRUE,mu=mu,sig=sig,prior=prior)
  #define posterior distribution
  }else{logpost=function(start,data) posterior(start,data)}
  if(is.null(propob)){
    propob = lapl_aprx(data[,1],data[,-1])
  }
  varprop = vscale*propob$var
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
  AcceptanceRate = AccptRate/iter
  val = list(Matpram=Mat[-c(1:burn),],AcceptanceRate=AcceptanceRate)
  cat("Random Walk MH algorithm successful. Acceptance ratio = ", AcceptanceRate," \n")
  return(val)
  }


#===============================================================================================#
#' Independence Metropolis-Hastings Algorithm
#'
#' \code{IndepMH} computes random draws of parameters using a specified proposal distribution.
#'
#' @param data data required for the posterior distribution
#' @param propob a list of mean and variance-covariance of the normal proposal distribution; 
#' defaults to \code{NULL}. Write list as propob=list(mode=mode,var=variance-covariance)
#' @param posterior the posterior distribution. It is set to null in order to use the logit posterior.
#' The user can specify log posterior as a function of parameters and data (pars,data). For a more
#' flexible and generic implementation, use \code{\link{indepMHgen}}
#' @param iter number of random draws desired (default: 1500)
#' @param burn burn-in period for the MH algorithm (default: 500)
#' @param vscale a positive value to scale up or down the variance-covariance matrix in
#' the proposal distribution
#' @param start starting values of parameters for the MH algorithm.
#' It is automatically generated but the user can also specify.
#' @param prior the prior distribution (default: "Normal", alternative: "Uniform")
#' @param mu the mean of the normal prior distribution (default:0)
#' @param sig the variance of the normal prior distribution (default:10)
#' @return val a list of matrix of draws pardraws and the acceptance rate
#'
#' @examples
#' y = indicat(faithful$waiting,70)
#' x = scale(cbind(faithful$eruptions,faithful$eruptions^2))
#' data = data.frame(y,x); propob<- lapl_aprx(y,x)
#' IndepMH_n<- IndepMH(data=data,propob,iter = 102, burn = 2) # prior="Normal"
#' IndepMH_u<- IndepMH(data=data,propob,prior="Uniform",iter = 102, burn = 2) # prior="Uniform"
#' par(mfrow=c(3,1));invisible(apply(IndepMH_n$Matpram,2,function(x)plot(density(x))))
#' invisible(apply(IndepMH_u$Matpram,2,function(x)plot(density(x))));par(mfrow=c(1,1))
#'
#' @export

IndepMH<- function(data,propob=NULL,posterior=NULL,iter=1500,burn=500,vscale=1.5,
                   start=NULL,prior="Uniform",mu=0,sig=10){
  if(is.null(posterior)){
    logpost<- function(start,data) posterior(start,data,Log=T,mu=mu,sig=sig,prior=prior)
    #define posterior distribution
  }else{logpost=function(start,data) posterior(start,data)}
  if(is.null(propob)){
    propob = lapl_aprx(data[,1],data[,-1])
  }
  varprop = vscale*propob$var 
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
  Accept_Rate = AccptRate/iter
  val = list(Matpram=Mat[-c(1:burn),],Accept_Rate = AccptRate/iter)
  cat("IndepMH algorithm successful. Acceptance ratio = ", Accept_Rate," \n")
  return(val)
}
#===========================================================================================#

#===========================================================================================#
#' A Generic Independence Metropolis-Hastings Algorithm
#'
#' \code{IndepMHgen} computes random draws of parameters using a normal proposal distribution.
#' This function implements a generic form of \code{\link{IndepMH}}
#'
#' @param start starting values of parameters for the MH algorithm.
#' It is automatically generated from the normal proposal distribution but the user can also specify.
#' @param posterior the log posterior distribution function.
#' Should take parameter input of the same length as \code{start} or \code{propob$mode} 
#' @param ... additional arguments to the posterior function
#' @param propob a list of mode and variance-covariance matrix of the normal proposal distribution. 
#' Write list as propob=list(mode=mode,var=variance-covariance)
#' @param const a vector function of parameters showing non-negative inequality constraints to be satisfied. 
#' @param seed an integer as seed for reproducibility
#' @param vscale a value multiplied by \code{propob$var} in order to adjust the proposal distribution. Else
#' set to the character string "HS18" for the Herbst and Shorfheide (2018) scale updating for a 0.25
#' acceptance ratio.
#' The default is \code{1.5} but the user can adjust it until a satisfactory acceptance rate is obtained.
#' @param iter number of random draws desired (default: 5000)
#' @param burn burn-in period for the MH algorithm (default: floor(iter/10))
#' @param report a numeric frequency (i.e. after how many iterations to report progress) for reporting
#'  algorithm progress; default - NULL
#' @return Matpram a matrix of parameter draws
#' @return postvals vector of posterior values corresponding to parameter draws \code{Matpram}
#' @return AcceptRatio the acceptance ratio
#' 
#' @examples 
#' #a toy example for illustration
#' ## f(c) = 1/(3.618*sqrt(pi))* exp(-0.6*(c[1]-2)^2-0.4*(c[2]+2)^2) 
#' # an improper posterior
#' logpost = function(c) -0.6*(c[1]-2)^2-0.4*(c[2]+2)^2 #log posterior distribution
#' optp<-optim(par=c(0,0),fn=logpost,control=list(fnscale=-1),hessian = TRUE) 
#' # laplace approximation of the posterior
#' propob = list(mode=optp$par,var=-solve(optp$hessian)) #parameters of proposal distribution
#' eigen(propob$var)$values # var-cov of proposal distribution is positive definite
#' MHobj<- indepMHgen(posterior = logpost,propob = propob,vscale = "HS18",iter = 6000,report=1000)
#' # create an independent Metropolis-Hastings object
#' dim(MHobj$Matpram) # a 2 x 5000 matrix with columns corresponding to draws of c1 and c2
#' par(mfrow=c(1,2))
#' hist(MHobj$Matpram[1,],20,main = "Histogram c1",xlab = "c1")
#' hist(MHobj$Matpram[2,],20,main = "Histogram c2",xlab = "c2"); par(mfrow=c(1,1))
#' MHobj$AcceptRatio # acceptance ratio
#'
#' @export

indepMHgen<- function(start=NULL,posterior=NULL,...,propob=NULL,const=NULL,
                      seed=1,vscale=1.5,iter=5000,burn=floor(0.1*iter),
                      report=NULL){
  if(vscale!="HS18"){
    varprop = vscale*propob$var
  }else{
    varprop = propob$var
    c0 = 1.0 #initialise adapting scale
  }
  if(!is.null(seed)){set.seed(seed = seed)}
  npar = length(propob$mode)
  Mat = array(0, c(iter, npar)); postvals<- c(0)
  
  if(is.null(const))
  {  
    if(is.null(start)){
      start = MASS::mvrnorm(n=1,propob$mode,varprop)
    }
    Mat[1,] = start; AccptRate<-0; postvals[1]<- posterior(start,...)
    
    for(i in 2:iter){
      start= Mat[i-1,]
      prop = MASS::mvrnorm(n=1,propob$mode,varprop)#make a draw from proposal dist
      lpa = posterior(prop,...); lpb = postvals[i-1]
      accprob = exp(lpa-lpb)
      if(is.na(accprob)){accprob=0} #penalise NA fun values
      # the other part cancels out because the normal distribution is symmetric
      if(stats::runif(1)< accprob){
        Mat[i,]=prop
        AccptRate<- AccptRate +1
        postvals[i]<- lpa
      }else{
        Mat[i,]=start 
        postvals[i]<- postvals[i-1]
      }
      if(vscale=="HS18"){
        rx = AccptRate/i #acceptance ratio so far
        c1 = c0*(0.95 + 0.1*exp(16*(rx-0.25))/(1+exp(16*(rx-0.25))))
        c0 = c1
      }
      if(!is.null(report)){
        if(i%%report==0){message(i," iterations done. ",I(iter-i)," more to go. \n")}
      }
    }
  }else{
    if(is.null(start)){
      cnt<- 0
      start = MASS::mvrnorm(n=1,propob$mode,varprop)
      while(any(const(start)<0)){
        start = MASS::mvrnorm(n=1,propob$mode,varprop) ; cnt<- cnt+1
        if(cnt>burn){stop("Cannot generate proposal draws within the constraint region")}
      }
    }
    Mat[1,] = start; AccptRate<-0; postvals[1]<- posterior(start,...)
    
    for(i in 2:iter){
      start= Mat[i-1,]
      prop = MASS::mvrnorm(n=1,propob$mode,varprop)#make a draw from proposal dist
      cnt<- 0
      while(any(const(prop)<0)){
        prop = MASS::mvrnorm(n=1,propob$mode,varprop) ; cnt<- cnt+1
        if(cnt>burn){stop("Cannot generate proposal draws within the constraint region")}
      }
      lpa = posterior(prop,...); lpb = postvals[i-1]
      accprob = exp(lpa-lpb)
      # the other part cancels out because the normal distribution is symmetric
      if(is.na(accprob)){accprob=0} #penalise NA fun values
      if(stats::runif(1)< accprob){
        Mat[i,]=prop
        AccptRate<- AccptRate +1
        postvals[i]<- lpa
      }else{
        Mat[i,]=start 
        postvals[i]<- postvals[i-1]
      }
      if(vscale=="HS18"){
        rx = AccptRate/i #acceptance ratio so far
        c1 = c0*(0.95 + 0.1*exp(16*(rx-0.25))/(1+exp(16*(rx-0.25))))
        c0 = c1
      }
      if(!is.null(report)){
        if(i%%report==0){message(i," iterations done. ",I(iter-i)," more to go. \n")}
      }
    }
  }
  if(report){message("indepMHgen algorithm successful\n")}
  val = list(Matpram=t(Mat[-c(1:burn),]),postvals=postvals[-c(1:burn)],AcceptRatio = AccptRate/iter)
  return(val)
}
#===========================================================================================#


#===========================================================================================#
#' A Generic Random Walk Metropolis-Hastings Algorithm
#'
#' \code{rwMHgen} computes random draws of parameters using a normal proposal distribution.
#' This function implements a generic form of \code{\link{RWMH}}
#'
#' @param start starting values of parameters for the MH algorithm.
#' It is automatically generated from the normal proposal distribution but the user can also specify.
#' @param posterior the log posterior distribution function. 
#' Should take parameter input of the same length as \code{start} or \code{propob$mode}
#' @param ... additional arguments to the posterior function
#' @param propob a list of mode and variance-covariance matrix of the normal proposal distribution. 
#' Write list as propob=list(mode=mode,var=variance-covariance)
#' @param const a vector function of parameters showing non-negative inequality constraints to be satisfied. 
#' @param seed an integer as seed for reproducibility
#' @param vscale a value multiplied by \code{propob$var} in order to adjust the proposal distribution. Else
#' set to the character string "HS18" for the Herbst and Shorfheide (2018) scale updating for a 0.25
#' acceptance ratio.
#' The default is \code{1.5} but the user can adjust it until a satisfactory acceptance rate is obtained.
#' @param iter number of random draws desired (default: 5000)
#' @param burn burn-in period for the MH algorithm (default: floor(0.1*iter))
#' @param report a numeric frequency (i.e. after how many iterations to report progress) for reporting
#'  algorithm progress; default - NULL
#' @return Matpram a matrix of parameter draws
#' @return postvals vector of posterior values corresponding to parameter draws \code{Matpram}
#' @return AcceptRatio the acceptance ratio
#' 
#' @examples 
#' #a toy example for illustration
#' ## f(c) = 1/(3.618*sqrt(pi))* exp(-0.6*(c[1]-2)^2-0.4*(c[2]+2)^2) 
#' # an improper posterior
#' logpost = function(c) -0.6*(c[1]-2)^2-0.4*(c[2]+2)^2 #log posterior distribution
#' optp<-optim(par=c(0,0),fn=logpost,control=list(fnscale=-1),hessian = TRUE) 
#' # laplace approximation of the posterior
#' propob = list(mode=optp$par,var=-solve(optp$hessian)) #parameters of proposal distribution
#' eigen(propob$var)$values # var-cov of proposal distribution is positive definite
#' MHobj<- rwMHgen(posterior = logpost,propob = propob,vscale = "HS18",iter = 6000,report=1500)
#' # create an independent Metropolis-Hastings object
#' dim(MHobj$Matpram) # a 2 x 5000 matrix with columns corresponding to draws of c1 and c2
#' par(mfrow=c(1,2))
#' hist(MHobj$Matpram[1,],20,main = "Histogram c1",xlab = "c1")
#' hist(MHobj$Matpram[2,],20,main = "Histogram c2",xlab = "c2"); par(mfrow=c(1,2))
#' MHobj$AcceptRatio # acceptance ratio
#'
#' @export

rwMHgen<- function(start=NULL,posterior=NULL,...,propob=NULL,const=NULL,
                   seed=1,vscale=1.5,iter=5000,burn=floor(0.1*iter),
                   report=NULL){
  if(vscale!="HS18"){
    varprop = vscale*propob$var
  }else{
    varprop = propob$var
    c0 = 1.0 #initialise adapting scale
  }
  if(!is.null(seed)){set.seed(seed = seed)}
  npar = length(propob$mode)
  Mat = array(0, c(iter, npar)); postvals<- c(0)
  
  if(is.null(const))
  {  
    if(is.null(start)){
      start = MASS::mvrnorm(n=1,propob$mode,varprop)
    }
    Mat[1,] = start; AccptRate<-0; postvals[1]<- posterior(start,...)
    
    for(i in 2:iter){
      start= Mat[i-1,]
      prop = MASS::mvrnorm(n=1,start,varprop)#make a draw from proposal dist
      lpa = posterior(prop,...); lpb = postvals[i-1]
      accprob = exp(lpa-lpb)
      if(is.na(accprob)){accprob=0} #penalise NA fun values
      # the other part cancels out because the normal distribution is symmetric
      if(stats::runif(1)< accprob){
        Mat[i,]=prop
        AccptRate<- AccptRate +1
        postvals[i]<- lpa
      }else{
        Mat[i,]=start 
        postvals[i]<- postvals[i-1]
      }
      if(vscale=="HS18"){
        rx = AccptRate/i #acceptance ratio so far
        c1 = c0*(0.95 + 0.1*exp(16*(rx-0.25))/(1+exp(16*(rx-0.25))))
        c0 = c1
      }
      if(!is.null(report)){
        if(i%%report==0){message(i," iterations done. ",I(iter-i)," more to go. \n")}
      }
    }
  }else{
    if(is.null(start)){
      cnt<- 0
      start = MASS::mvrnorm(n=1,propob$mode,varprop)
      while(any(const(start)<0)){
        start = MASS::mvrnorm(n=1,propob$mode,varprop) ; cnt<- cnt+1
        if(cnt>burn){stop("Cannot generate proposal draws within the constraint region")}
      }
    }
    Mat[1,] = start; AccptRate<-0; postvals[1]<- posterior(start,...)
    
    for(i in 2:iter){
      start= Mat[i-1,]
      prop = MASS::mvrnorm(n=1,start,varprop)#make a draw from proposal dist
      cnt<- 0
      while(any(const(prop)<0)){
        prop = MASS::mvrnorm(n=1,start,varprop) ; cnt<- cnt+1
        if(cnt>burn){stop("Cannot generate proposal draws within the constraint region")}
      }
      lpa = posterior(prop,...); lpb = postvals[i-1]
      accprob = exp(lpa-lpb)
      # the other part cancels out because the normal distribution is symmetric
      if(is.na(accprob)){accprob=0} #penalise NA fun values
      if(stats::runif(1)< accprob){
        Mat[i,]=prop
        AccptRate<- AccptRate +1
        postvals[i]<- lpa
      }else{
        Mat[i,]=start 
        postvals[i]<- postvals[i-1]
      }
      if(vscale=="HS18"){
        rx = AccptRate/i #acceptance ratio so far
        c1 = c0*(0.95 + 0.1*exp(16*(rx-0.25))/(1+exp(16*(rx-0.25))))
        c0 = c1
      }
      if(!is.null(report)){
        if(i%%report==0){message(i," iterations done. ",I(iter-i)," more to go. \n")}
      }
    }
  }
  if(report){message("rwMHgen algorithm successful\n")}
  val = list(Matpram=t(Mat[-c(1:burn),]),postvals=postvals[-c(1:burn)],AcceptRatio = AccptRate/iter)
  return(val)
}
#===========================================================================================#