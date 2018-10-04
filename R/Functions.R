#' Indicator function
#'
#' This function creates 0-1 indicators for a given threshold y0 and vector y
#'
#' @param y vector y
#' @param y0 threshold value y0
#' @return val
#'
#' @export

indicat<- function(y,y0){
  if(y<=y0){
    val<- 1
  }else{
    val=0
  }
  return(val)
}
indicat<- Vectorize(indicat) #make the function usable with a vector y

#===============================================================================
#' Logit function
#'
#' This is the link function for logit regression
#'
#' @param x Random variable
#' @return val Probability value from the logistic function
#'
#' @export

LogitLink = function(x) {
  val=(1/(1+exp(-x)))
  return(val)
  }
LogitLink=Vectorize(LogitLink) # ensure usability with entire vectors

#================================================================================================#
#' Logit likelihood function
#'
#' \code{logit} is the logistic likelihood function given data.
#'
#' @param start starting values
#' @param data dataframe. The first column should be the dependent variable.
#' @param Log a logical input (default True) to take the log of the likelihood. Only set 
#' FALSE if you do not need the default TRUE
#' @return like returns the likelihood function value.
#'
#' @examples
#' # vv<- logit(start=rep(0,5),data = genmle::dat_mroz) # log-likelihood
#' # vv1<- logit(start=rep(0,5),data = genmle::dat_mroz,Log = F)# likelihood
#'
#' @export

logit = function(start,data,Log=TRUE){
  y=data[,1]; x=data[,-1]; #x = scale(x) # apportion response and explanatory variables
  penalty <- -10^12
  like <- penalty
  N <- length(y);
  npar <- length(start)
  f <- matrix(0,N,1)
  x = as.matrix(data.frame(rep(1,N),x)) # add ones for intercept term
  m = c(0);
  for (i in 1:nrow(data)) {
    m[i] = x[i,]%*%start
  }
  g = LogitLink(m) # take logit transform of all values
  # multiply g by weights if weights
  for (i in 1:nrow(data)) {
    if(y[i]==1){p=g[i]}
    if(y[i]==0){p=(1-g[i])}
    f[i] = ifelse(Log, log(p), p)  #compute log of likelihood
  }
  like = ifelse(Log, sum(f),prod(f))
  if (is.na(like) | is.infinite(like)) {
    like <- penalty
  }
  return(like) # for maximisation
}
#================================================================================================#
#' Normal Prior distribution
#'
#' This normal prior distribution is a product of univariate N(mu,sig)
#'
#' @param pars parameter values
#' @param mu mean value of each parameter value
#' @param sig standard deviation of each parameter value
#' @param Log logical to take the log of likelihoods or not (default: FALSE)
#' @return val Product of probability values for each parameter
#'
#' @examples
#' # prior_n(rep(0,6),0,10,Log = T) #log of prior
#' # prior_n(rep(0,6),0,10,Log = F) #no log
#'
#' @export
prior_n<- function(pars,mu,sig,Log=F){
  val<- ifelse(!Log, prod(stats::dnorm(pars,mu,sig)),sum(stats::dnorm(pars,mu,sig,log= T))) 
  #if log FALSE
  return(val)
  }

#================================================================================================#
#' Uniform Prior distribution
#'
#' This uniform prior distribution proportional to 1 (This function is redundant for now)
#'
#' @param pars parameter values
#' @return val value of joint prior =1 for the uniform prior
#'
#' @export
prior_u<- function(pars){
  val<- 1
  return(val)
}


#================================================================================================#
#' Posterior distribution
#'
#' \code{posterior} computes the value of the posterior at parameter values pars
#'
#' @param pars parameter values
#' @param data dataframe. The first column must be the binary dependent variable
#' @param Log logical to take the log or not (default: TRUE)
#' @param mu mean of prior of each parameter value in case the prior is Normal (default: 0)
#' @param sig standard deviation of prior of each parameter in case the prior is Normal 
#' (default: 25)
#' @param prior string input of "Normal" or "Uniform" prior distribution to use
#' @return val value function of the posterior
#'
#' @examples
#' # posterior(rep(0,5),genmle::dat_mroz,Log = F,mu=0,sig = 10,prior = "Normal") # no log
#' # posterior(rep(0,5),genmle::dat_mroz,Log = T,mu=0,sig = 10,prior = "Normal") # log
#' # posterior(rep(0,5),genmle::dat_mroz,Log = T) # use default values
#'
#' @export
posterior<- function(pars,data,Log=T,mu=0,sig=25,prior="Normal"){
  if(Log){
    val = logit(pars,data,Log=T) +
      ifelse(identical("Normal",prior),prior_n(pars,mu,sig,Log=T),0)
  }else{
    val = logit(pars,data,Log=F) *
      ifelse(identical("Normal",prior),prior_n(pars,mu,sig,Log=F),1)
  }

  return(val)
}

#================================================================================================#
#' Laplace approximation of posterior to normal
#'
#' This function generates mode and variance-covariance for a normal proposal
#' distribution for the bayesian logit.
#'
#' @param y the binary dependent variable y
#' @param x the matrix of independent variables.
#' @param glmobj logical for returning the logit glm object
#' @return val A list of mode variance-covariance matrix, and scale factor for
#' proposal draws from the multivariate normal distribution.
#'
#' @examples
#' # dat=genmle::dat_mroz ;y=dat$y; x = dat[,-1]; gg<- lapl_aprx(y,x)
#'
#' @export
lapl_aprx<- function(y,x,glmobj=FALSE){ #laplace approximation
  dat = data.frame(y,x)
  lgitob<-stats::glm(dat$y~.,data=dat,family = "binomial")
  
  if(!glmobj){
  val<- list(mode= lgitob$coefficients,var = stats::vcov(lgitob))
  }else{
    val<- list(mode= lgitob$coefficients,var = stats::vcov(lgitob)
               ,glmobj=lgitob)
  }
  return(val)
}
#================================================================================================#
#' An \code{optim}-based function maximiser
#'
#' \code{lapl_aprx2} is used to take a laplace approximation of the posterior distribution
#'
#' @param start the starting parameter values
#' @param func the function to be maximised
#' @param ... other arguments to be passed to func
#' @return val list of mode and variance-covariance matrix
#'
#' @export

lapl_aprx2<- function(start,func,...){
  mxob<- stats::optim(par=start,fn=func,...,
               control=list(fnscale=-1,trace=F,REPORT=50,
                            factr=1e-15, maxit=20000,
                            ndeps=rep(1e-6,length(start))),hessian=T)
  val = list(mode=mxob$par,var=solve(-mxob$hessian))
  return(val)
}


#================================================================================================#
#' Fitted logit probabilities
#'
#' \code{fitlogit} obtains a vector of fitted logit probabilities given parameters (pars)
#' and data
#'
#' @param pars vector of parameters
#' @param data data frame. The first column of the data frame ought to be the binary dependent 
#' variable
#' @return vec vector of fitted logit probabilities
#'
#' @export


fitlogit<- function(pars,data){
  y=data[,1]; x=data[,-1]; #x = scale(x) # apportion response and explanatory variables
  penalty <- -10^12
  like <- penalty
  N <- length(y);
  npar <- length(pars)
  f <- matrix(0,N,1)
  x = as.matrix(data.frame(rep(1,N),x)) # add ones for intercept term
  m = c(0);
  for (i in 1:nrow(data)) {
    m[i] = x[i,]%*%pars
  }
  vec = LogitLink(m) # take logit transform of all values
  return(vec)
}

#================================================================================================#
#' The distribution of mean fitted logit probabilities
#'
#' \code{fitdist} function generates a vector of mean fitted probabilities that constitute the 
#' distribution
#'
#' @param Matparam an M x k matrix of parameter draws, each being 1 x k
#' @param data dataframe used to obtain Matparam
#' @return dist fitted distribution
#'
#' @export
fitdist<- function(Matparam,data){
  fm = function(i) mean(fitlogit(Matparam[i,],data))
  dist<-sapply(1:nrow(Matparam),fm)
  return(dist)
}

#================================================================================================#
#' Parlapply a function
#'
#' \code{paRLply} parlapply from the \code{parallel} package with a function as input
#'
#' @param vec vector of inputs over which to parallel compute
#' @param fn the function
#' @param type this option is set to "FORK", use "PSOCK" for windows
#' @param ... extra inputs to \code{fn}
#' @return out parallel computed output
#'
#' @export
parLply<- function(vec,fn,type="FORK",...){
no_cores<-parallel::detectCores() - 1
c1<-parallel::makeCluster(no_cores, type = type)
out<- parallel::parLapply(c1,vec,fn,...)
parallel::stopCluster(c1)
out
}
#================================================================================================#
#' Quantile conversion of a bayesian distribution matrix
#'
#' \code{quant_bdr} converts a bayesian distribution regression matrix from \code{par_distreg()}
#' output to a matrix of quantile distribution.
#'
#' @param taus a vector of quantile indices
#' @param thresh a vector of threshold values used in a \code{par_distreg()} type function
#' @param mat bayesian distribution regression output matrix
#' @return qmat matrix of quantile distribution
#'
#' @export
quant_bdr<- function(taus,thresh,mat){
  mat=t(apply(mat,1,sort)); 
  q.inv<-function(j){
    zg=stats::spline(x=sort(mat[,j]),y=thresh,xout = taus)
    sort(zg$y)
  } ; q.inv=Vectorize(q.inv)
  qmat<-sapply(1:ncol(mat),q.inv)
  qmat
}

#=====================================================================================================#
#' Smooth marginal conditional distribution
#'
#' \code{getsmooth} uses moving means to minimise kinks in marginal conditional distributions
#'
#' @param vec a vector of representing the marginal conditional quantile
#' @param Step length of moving window used to calculate the moving mean
#' @return val vector of smoothed marginal conditional quantile
#'
#' @export

getsmooth<- function(vec,Step=NULL){
  #vec<- sort(vec); 
  lv<- length(vec)
  d<- rep(0,lv)
  if(is.null(Step)){Step=3}
  d[1:Step]<- cumsum(vec[1:Step])/c(1:Step)
  if(lv>Step){
    for(j in (Step+1):lv){
      d[j]<- mean(vec[(j-Step+1):j])
    }
  }else{
    stop("The length of the vector cannot be longer than Step")
  }
  #val<- sort(d)
  val<- d
  return(val)
}


#=====================================================================================================#
#' Smooth marginal conditional distribution using smooth.spline()
#'
#' \code{getsmooth.spline} uses the function \code{smooth.spline()} in the 
#' \code{stats} package for smoothing
#'
#' @param x a vector of representing the marginal conditional quantile
#' @param y length of moving window used to calculate the moving mean
#' @param CI (optional) a list of \code{CI$lb} and \code{CI$up} of lower and upper confidence
#' intervals of \code{y}
#' @param cv logical for cross-validation to be passed to \code{smooth.spline()}
#' @param ... additional paramters to be passed to \code{smooth.spline()}
#' @return out a list with the following components:
#' @return x output vector x
#' @return y output vector y
#' @return lb lower confidence band if \code{CI} is provided
#' @return ub upper confidence band if \code{CI} is provided
#'
#' @export

getsmooth.spline<- function(x,y,CI=NULL,cv=TRUE,...){
  hz<- stats::smooth.spline(x=x,y=y,cv=cv,...)
  if(is.null(CI)){
  out=list(x=hz$x,y=hz$y)
  }else{
    # hl<- smooth.spline(x=x,y=CI$lb,cv=cv,...)
    # hu<- smooth.spline(x=x,y=CI$ub,cv=cv,...)
    up<- mean(abs(CI$ub-y)); lw<- mean(abs(CI$lb-y))
    lb= hz$y-lw; ub=hz$y+up
    #out=list(x=hz$x,y=hz$y,lb=hl$y,ub=hu$y)
    out=list(x=hz$x,y=hz$y,lb=lb,ub=ub)
  }
}

#=====================================================================================================#
#' Simultaneous bayesian confidence bands
#'
#' \code{simcnfB} obtains bayesian distribution confidence bands
#'
#' @param DF the target distribution as a vector
#' @param DFmat the matrix of draws of the distribution, rows correspond to 
#' indices elements in \code{DF}
#' @param alpha level such that \code{1-alpha} is the desired probability of coverage
#' @return cstar - a constant to add and subtract from DF to create 
#' the confidence bands.
#' @return uCB - upper band
#'
#' @export
#' 
simcnfB<- function(DF,DFmat,alpha=0.05){
  M = abs(DFmat-DF)
  #sdM = apply(M, 1, sd)
  sdM = apply(M, 1, quantile,probs=0.75)-apply(M, 1, quantile,probs=0.25)
  Ms = M/sdM
  cstar<-max(apply(Ms,1,quantile,probs=(1-alpha)))
  #cstar<-quantile(apply(Ms,1,max),probs=(1-alpha))
  cstar*sdM
}
















