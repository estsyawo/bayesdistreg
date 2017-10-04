#' An optim-based function maximiser
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
  mxob<- optim(par=start,fn=func,...,
               control=list(fnscale=-1,trace=F,REPORT=50,
                            factr=1e-15, maxit=20000,
                            ndeps=rep(1e-6,length(start))),hessian=T)
  val = list(mode=mxob$par,var=solve(-mxob$hessian))
  return(val)
  }
