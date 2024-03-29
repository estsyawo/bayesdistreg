#' Indicator function
#'
#' This function creates 0-1 indicators for a given threshold y0 and vector y
#'
#' @param y vector y
#' @param y0 threshold value y0
#' @return val
#' 
#' @examples
#' yb = indicat(faithful$waiting,quantile(faithful$waiting,probs=c(0.25,0.5,0.75))) 
#' summary(yb)
#' @export

indicat<- function(y,y0){
  fn=function(y0){as.numeric(y<=y0)}
  fn=Vectorize(fn)
  sapply(y0,fn)
}
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
LogitLink=base::Vectorize(LogitLink) # ensure usability with entire vectors

#================================================================================================#
#' Logit likelihood function
#'
#' \code{logit} is the logistic likelihood function given data.
#'
#' @param start vector of starting values
#' @param data dataframe. The first column should be the dependent variable.
#' @param Log a logical input (defaults to \code{True}) to take the log of the likelihood. 
#' @return like returns the likelihood function value.
#' 
#' @examples
#' y = indicat(faithful$waiting,mean(faithful$waiting)) 
#' x = scale(cbind(faithful$eruptions,faithful$eruptions^2))
#' data = data.frame(y,x)
#' logit(rep(0,3),data)
#' @export

logit = function(start,data,Log=TRUE){
  y=data[,1]; x=data[,-1]; # apportion response and explanatory variables
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
  g = LogitLink(m) # take logit transform of linear prediction vector
  for (i in 1:nrow(data)) {
    if(y[i]==1){p=g[i]}
    if(y[i]==0){p=(1-g[i])}
    f[i] = ifelse(Log, log(p), p)  # compute log of likelihood
  }
  like = ifelse(Log, sum(f),prod(f))
  if (is.na(like) | is.infinite(like)) {
    like <- penalty
  }
  return(like)
}
#================================================================================================#
#' Normal Prior distribution
#'
#' This normal prior distribution is a product of univariate N(mu,sig)
#'
#' @param pars parameter values
#' @param mu mean value of each parameter value
#' @param sig standard deviation of each parameter value
#' @param Log logical to take the log of prior or not (defaults to FALSE)
#' @return val Product of probability values for each parameter
#'
#' @examples
#' prior_n(rep(0,6),0,10,Log = TRUE) #log of prior
#' prior_n(rep(0,6),0,10,Log = FALSE) #no log
#'
#' @export
prior_n<- function(pars,mu,sig,Log=FALSE){
  val<- ifelse(!Log, prod(stats::dnorm(pars,mu,sig)),sum(stats::dnorm(pars,mu,sig,log= T))) 
  #if log FALSE
  return(val)
  }

#================================================================================================#
#' Uniform Prior distribution
#'
#' This uniform prior distribution proportional to 1
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
#' \code{posterior} computes the value of the posterior at parameter values \code{pars}
#'
#' @param pars parameter values
#' @param data dataframe. The first column must be the binary dependent variable
#' @param Log logical to take the log of the posterior.(defaults to TRUE)
#' @param mu mean of prior of each parameter value in case the prior is Normal (default: 0)
#' @param sig standard deviation of prior of each parameter in case the prior is Normal 
#' (default: 25)
#' @param prior string input of "Normal" or "Uniform" prior distribution to use
#' @return val value function of the posterior
#'
#' @examples
#' y = indicat(faithful$waiting,mean(faithful$waiting)) 
#' x = scale(cbind(faithful$eruptions,faithful$eruptions^2))
#' data = data.frame(y,x)
#' posterior(rep(0,3),data,Log = FALSE,mu=0,sig = 10,prior = "Normal") # no log
#' posterior(rep(0,3),data,Log = TRUE,mu=0,sig = 10,prior = "Normal") # log
#' posterior(rep(0,3),data,Log = TRUE) # use default values
#'
#' @export
posterior<- function(pars,data,Log=TRUE,mu=0,sig=25,prior="Normal"){
  if(Log){
    val = logit(pars,data,Log=TRUE) +
      ifelse(identical("Normal",prior),prior_n(pars,mu,sig,Log=T),0)
  }else{
    val = logit(pars,data,Log=FALSE) *
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
#'  y = indicat(faithful$waiting,mean(faithful$waiting)) 
#'  x = scale(cbind(faithful$eruptions,faithful$eruptions^2))
#'  gg<- lapl_aprx(y,x)
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
#' Laplace approximation of posterior to normal
#'
#' \code{lapl_aprx2} is a more flexible alternative to \code{lapl_aprx}. This creates
#' \code{glm} objects from which joint asymptotic distributions can be computed.
#'
#' @param y the binary dependent variable y
#' @param x the matrix of independent variables.
#' @param family a parameter to be passed \code{\link[stats]{glm}}, defaults to the logit model
#' @param ... additional parameters to be passed to \code{\link[stats]{glm}}
#' @return val A list of mode variance-covariance matrix, and scale factor for
#' proposal draws from the multivariate normal distribution.
#' 
#' @examples
#' y = indicat(faithful$waiting,mean(faithful$waiting)) 
#' x = scale(cbind(faithful$eruptions,faithful$eruptions^2))
#' (gg<- lapl_aprx2(y,x)); coef(gg); vcov(gg)
#'
#' @export

lapl_aprx2<- function(y,x,family = "binomial",...){ #laplace approximation
  dat = data.frame(y,x)
  lgitob<-stats::glm(dat$y~.,data=dat,family = family,...)
  lgitob
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
  y=data[,1]; x=data[,-1]; 
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
#' distribution. This involves marginalising out covariates.
#'
#' @param Matparam an M x k matrix of parameter draws, each being a 1 x k vector
#' @param data dataframe used to obtain Matparam
#' @return dist fitted (marginalised) distribution
#'
#' @export
fitdist<- function(Matparam,data){
  fm = function(i) mean(fitlogit(Matparam[i,],data))
  dist<-sapply(1:nrow(Matparam),fm)
  return(dist)
}

#================================================================================================#
#' Parallel compute
#'
#' \code{parLply} uses \code{\link[parallel]{parLapply}} with a function as input
#'
#' @param vec vector of inputs over which to parallel compute
#' @param fn the function
#' @param type this option is set to "FORK", use "PSOCK" on windows
#' @param no_cores the number of cores to use. Defaults at 1
#' @param ... extra inputs to \code{fn()}
#' @return out parallel computed output
#'
#' @export
parLply<- function(vec,fn,type="FORK",no_cores=1,...){
c1<-parallel::makeCluster(no_cores, type = type)
out<- parallel::parLapply(c1,vec,fn,...)
parallel::stopCluster(c1)
out
}
#================================================================================================#
#' Quantile conversion of a bayesian distribution matrix
#'
#' \code{quant_bdr} converts a bayesian distribution regression matrix from \code{\link{par_distreg}}
#' output to a matrix of quantile distribution.
#'
#' @param taus a vector of quantile indices
#' @param thresh a vector of threshold values used in a \code{\link{par_distreg}} type function
#' @param mat bayesian distribution regression output matrix; columns should correspond to draws of
#' the distribution function
#' @return qmat matrix of quantile distribution
#' 
#' @examples 
#' set.seed(2); Fmat = matrix(runif(1000),1000,10); taus = seq(0,1, length.out=1000)
#' thresh = qnorm(taus)
#' par(mfrow=c(1,2))
#' plot(thresh,sort(apply(apply(Fmat,2,sort),1,mean)),type="l",xlab="y",ylab="F(y)",
#' main="Distribution Function")
#' qmat<- quant_bdr(taus,thresh,Fmat); plot(taus,apply(qmat,1,mean),type = "l",
#' xlab="tau",ylab="Qy",main="Quantile Function")
#' par(mfrow=c(1,1))
#' @export

quant_bdr<- function(taus,thresh,mat){
  # define function internally
  f2q=function(fs,y0,qs){
    qss=rep(0,length(qs))
    for (i in 1:length(qs)) {
      fs=fs[!is.na(fs)] # remove NA in a vector
      if(max(fs) >= qs[i]){
        qss[i]= y0[which(sort(fs) >= qs[i])[1]]
      }
      else {qss[i]=max(y0)}
    }
    qss
  }
  # define function internally
  ffq<- function(j) f2q(fs=mat[,j],y0=thresh,qs=taus)
  qmato<- sapply(1:ncol(mat), ffq)
  qmato
}

#=====================================================================================================#
#' Symmetric simultaneous bayesian confidence bands
#'
#' \code{simcnfB} obtains symmetric bayesian distribution confidence bands
#'
#' @param DF the target distribution/quantile function as a vector
#' @param DFmat the matrix of draws of the distribution, rows correspond to 
#' elements in \code{DF}
#' @param alpha level such that \code{1-alpha} is the desired probability of coverage
#' @param scale logical for scaling using the inter-quartile range
#' @return cstar - a constant to add and subtract from DF to create 
#'  confidence bands if no scaling=FALSE else a vector of length DF.
#'  
#' @examples 
#' set.seed(14); m=matrix(rbeta(500,1,4),nrow = 5) + 1:5
#' DF = apply(m,1,mean); plot(1:5,DF,type="l",ylim = c(0,max(m)), xlab = "Index")
#' symCB<- simcnfB(DF,DFmat = m)
#' lines(1:5,DF-symCB,lty=2); lines(1:5,DF+symCB,lty=2)
#'
#' @export
#' 
simcnfB<- function(DF,DFmat,alpha=0.05,scale=FALSE){
  Ms = abs(DFmat-DF)
  if(!scale){
    dj<- apply(Ms,2,max,na.rm=TRUE)
    cstar<- stats::quantile(dj[!is.infinite(dj)],probs = (1-alpha),na.rm = TRUE)
    ans = cstar
  }else{
    scv = apply(Ms,1,stats::IQR,na.rm=TRUE)
    Msc = Ms/scv
    dj<- apply(Msc,2,max,na.rm=TRUE)
    cstar<- stats::quantile(dj[!is.infinite(dj)],probs = (1-alpha),na.rm = TRUE)
    ans = cstar*scv
  }
  
  return(ans)
}
#=====================================================================================================#
#' Asymmetric simultaneous bayesian confidence bands
#'
#' \code{asymcnfB} obtains asymmetric bayesian confidence bands on a target function
#'
#' @param DF the target distribution/quantile function as a vector
#' @param DFmat the matrix of draws of the distribution, rows correspond to 
#' elements in \code{DF}
#' @param alpha level such that \code{1-alpha} is the desired probability of coverage
#' @param scale logical for scaling using the inter-quartile range
#' @return a list containing 
#' \itemize{
#' \item cmin - constant to subtract from \code{DF} for the lower confidence band
#' \item cmax - constant to add to \code{DF} for the upper confidence band
#' }
#'  
#' @examples 
#' set.seed(14); m=matrix(rbeta(500,1,4),nrow = 5) + 1:5
#' DF = apply(m,1,mean); plot(1:5,DF,type="l",ylim = c(min(m),max(m)), xlab = "Index")
#' asyCB<- asymcnfB(DF,DFmat = m)
#' lines(1:5,DF-asyCB$cmin,lty=2); lines(1:5,DF+asyCB$cmax,lty=2)
#'
#' @export
#' 
asymcnfB<- function(DF,DFmat,alpha=0.05,scale=FALSE){
  alf2 = alpha/2
  Ms = DFmat-DF
  if(!scale){
    djmx<- apply(Ms,2,max,na.rm=TRUE)
    djmn<- apply(Ms,2,min,na.rm=TRUE)
    cmin = -stats::quantile(djmn,probs=alf2,na.rm = TRUE)
    cmax = stats::quantile(djmx,probs=(1-alf2),na.rm = TRUE)
    ans=list(cmin=cmin,cmax=cmax)
  }else{
    scv = apply(Ms,1,stats::IQR,na.rm=TRUE)
    Msc = Ms/scv
    djmx<- apply(Msc,2,max,na.rm=TRUE)
    djmn<- apply(Msc,2,min,na.rm=TRUE)
    cmin = -stats::quantile(djmn,probs=alf2,na.rm = TRUE)
    cmax = stats::quantile(djmx,probs=(1-alf2),na.rm = TRUE)
    ans=list(cmin=cmin*scv,cmax=cmax*scv)
  }
  return(ans)        
}
#===================================================================================================
#'  Montiel Olea and Plagborg-Moller (2018) confidence bands
#'
#' \code{jntCBOM} implements calibrated symmetric confidence bands (algorithm 2)
#' in Montiel Olea and Plagborg-Moller (2018).
#'
#' @param DF the target distribution/quantile function as a vector of length G
#' @param DFmat the matrix of draws of the distribution, rows correspond to 
#' elements in \code{DF}
#' @param alpha level such that \code{1-alpha} is the desired probability of coverage
#' @param eps steps of the grid on the interval (alpha/(2*G),alpha/2) to search
#' @return CB - confidence band, zeta - the optimal level
#' 
#' @examples 
#' set.seed(14); m=matrix(rbeta(500,1,4),nrow = 5) + 1:5
#' DF = apply(m,1,mean); plot(1:5,DF,type="l",ylim = c(min(m),max(m)), xlab = "Index")
#' jOMCB<- jntCBOM(DF,DFmat = m)
#' lines(1:5,jOMCB$CB[,1],lty=2); lines(1:5,jOMCB$CB[,2],lty=2)
#' 
#' @export
#' 
jntCBOM<- function(DF,DFmat,alpha=0.05,eps=1e-3){
  G = length(DF); DF = sort(DF)
  zeta = (alpha/(2*G))
  cp = 1;
  while(cp>(1-alpha) & zeta <=(alpha/2)){
    vm=apply(DFmat,1,stats::quantile,probs=c(zeta,(1-zeta)))
    vm=apply(vm,1,sort)
    zk = vm - DF
    cp=length(which(zk[,1]<=0 & zk[,2]>=0))/G
    zeta=zeta+eps
  }
  return(list(CB = vm,zeta=(zeta-eps)))
}

#===================================================================================================
#'  A Bayesian "p-value"
#'
#' \code{psimval} finds a minimum alpha in the interval (0,1) such that the 1-alpha simultaneous 
#' confidence band excludes the vector of DF. This is intended to be a summary value for a joint
#' hypothesis test of equality. 
#'
#' @param DF the target vector needed to construct the confidence band
#' @param DFmat the matrix of draws of the distribution, rows correspond to 
#' elements in \code{DF}
#' @param ialpha a starting value of alpha to be searched on the interval (0,1)
#' @param eps steps of the grid on the interval (0,1) to search
#' @param typeband the type of confidence band to use: "sym" for \code{\link{simcnfB}}, 
#' "asym" for \code{\link{asymcnfB}}, or "OM" for \code{\link{jntCBOM}} are supported.
#' @return alf the "p-value"
#' 
#' @examples 
#' set.seed(14); m=matrix(rbeta(500,1,4),nrow = 5) + 1:5
#' (pval<- round(psimval(DF=rep(0.1,5),DFmat = m,ialpha=0.05,typeband="OM"),3))
#' 
#' @export
#' 

psimval<-function(DF,DFmat,ialpha=0.05,eps=1e-04,typeband=c("sym","asym","OM")){
  typeband<- match.arg(typeband)
  if(typeband=="sym"){
    fn=function(DF,DFmat,alpha){
      c1a=simcnfB(DF=DF,DFmat = DFmat,alpha = alpha)  
      list(lbCB=(DF-c1a),ubCB=(DF+c1a))
    }
  }else if(typeband=="asym"){
    fn=function(DF,DFmat,alpha){
      c1a=asymcnfB(DF=DF,DFmat = DFmat,alpha = alpha)  
      list(lbCB=(DF-c1a$cmin),ubCB=(DF+c1a$cmax))
    }
  }else if(typeband=="OM"){
    fn=function(DF,DFmat,alpha){
      c1a=jntCBOM(DF=DF,DFmat = DFmat,alpha = alpha)  
      list(lbCB=(c1a$CB[,1]),ubCB=(c1a$CB[,2]))
    }
  }
  CB=fn(DF=DF,DFmat = DFmat,alpha = ialpha)
  
  zg = !all((CB$lbCB)*(CB$ubCB)<=0)
  alf = ialpha
  if(zg){ #deacrease alpha, else 
    while(zg & alf>eps){
      alf=alf-eps
      CB=fn(DF=DF,DFmat = DFmat,alpha = alf)
      zg = !all((CB$lbCB)*(CB$ubCB)<=0)
    }
    if(!zg){alf=alf+eps} #the last step violates; choose penultimate
  }else{#increase alpha
    while(!zg & alf<(1-eps)){
      alf=alf+eps
      CB=fn(DF=DF,DFmat = DFmat,alpha = alf)
      zg = !all((CB$lbCB)*(CB$ubCB)<=0)
    }
    if(zg){alf=alf-eps} #the last step violates; choose penultimate
  }
  return(alf)
}
