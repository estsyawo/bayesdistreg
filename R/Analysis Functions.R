#=====================================================================================================#
#' Compute bivariate distribution
#'
#' \code{fnCDF} generates the cumulative weights of a vector given a normalising constant
#'
#' @param vec vector of input
#' @param nrmc normalising constant
#' @return val a vector of normalised cumulative weights
#'
#' @examples
#' kd<- matrix(runif(12,0,10),ncol = 3); apply(kd,2,fnCDF,nrmv=max(apply(kd,2,sum)))
#' @export

fnCDF<- function(vec,nrmc=NULL){
  if(is.null(nrmc)){
    nrmc=sum(vec)
  }
  val<- cumsum(vec)/nrmc
  return(val)
}

#=====================================================================================================#
#' Generate weights for bayesian distribution regression 3D plot
#'
#' \code{cdfw} generates a matrix of weights required for dimension z of 3D plot
#'
#' @param mat matrix of draws from the bivariate distribution
#' @param lgyt the length of grid values of average fitted probabilities to use (default: 512)
#' @param nrmc the choice of normalising constant for the function \code{fnCDF} (default: NULL which corresponds to column sums of pdf weights, alternative: "MaxColSum" maximum column sum)
#' @return val a matix with the same dimensions as mat for the the 3D plot
#' @export
cdfw<- function(mat,lgyt= 512,nrmc=NULL){
  ny0 = nrow(mat); ndr<- ncol(mat)
  gyt<- seq(0,1,length.out = lgyt)
  mat1<- matrix(0,ny0,length(gyt))
  for(i in 1:ny0 ){
    gz<- t(mat)[,i]
    kz<- density(unname(gz))
    kzx=kz$x
    kzy=kz$y
    for(j in 1:length(gyt)){
      if(gyt[j]<min(kzx) | gyt[j]>max(kzx)){
        mat1[i,j]=0
      }else if(gyt[j]==min(kzx)){
        mat1[i,j]=kzy[1]
      }else{
        dd<-min(which(kzx>=gyt[j]))
        mat1[i,j]=mean(c(kzy[dd],kzy[dd-1]))
      }

    }
  }
  if(is.null(nrmc)){
    val=apply(mat1,2,fnCDF)
  }else{
    val=apply(mat1,2,fnCDF,nrmc=max(apply(mat1,2,sum)))
  }
  return(val)
}
#=====================================================================================================#
#' Get mode of a vector
#'
#' \code{getmode} is used to obtain the mode from a vector representing a constinuous distribution.
#' For a discrete variable, use the function \code{table()}
#'
#' @param vec vector of values from the continuous distribution
#' @return mval the mode of the continuous distribution
#'
#' @export
getmode<- function(vec){
  zz<- density(as.matrix(vec))
  mval<- zz$x[which.max(zz$y)]
  return(mval)
}


#=====================================================================================================#
#' Get marginal conditional distribution
#'
#' This function enables one to obtain the marginal conditional distribution at a
#' chosen distributional statistic, eg. mean, mode, median, or any quantile
#'
#' @param full_out full matrix of output. The rows correspond to threshold
#' values y_0 and the columns correspond to the number of MCMC draws
#' @param diststat distributional statistic; one can choose \code{"mean"}, \code{"median"}
#' \code{"mode"}, or \code{"quantile"}. The default is \code{NULL} which returns the median
#' @param tau the quantile. tau=0.5 corresponds to the median. this option is turned off
#' when \code{diststat} is set to another distributional statistic.
#' @return val this is a vector of the marginal conditional distribution
#'
#' @export
mcondist<- function(full_out,diststat=NULL,tau=0.5){
  full_out<-as.matrix(full_out)
  if(identical(diststat,"mean")){
    val<- apply(full_out,1,mean)
  }else if(identical(diststat,"median")){
    val<- apply(full_out,1,median)
  }else if(identical(diststat,"mode")){
    val<- apply(full_out,1,getmode)
  }else if(is.null(diststat) | identical(diststat,"quantile")){
    val<- apply(full_out,1,quantile,tau); val<- unname(val)
  }else{
    stop("Distributional statistic not reconised.")
  }
  return(sort(val))
}
#=====================================================================================================#
#' Smooth marginal conditional distribution
#'
#' \code{getsmooth} uses moving means to minimise kinks in marginal conditional distributions
#' at quantiles and modes
#'
#' @param vec a vector of representing the marginal conditional quantile
#' @param Step length of moving window used to calculate the moving mean
#' @return val vector of smoothed marginal conditional quantile
#'
#' @export

getsmooth<- function(vec,Step=NULL){
  vec<- sort(vec); lv<- length(vec)
  d<- rep(0,lv)
  if(is.null(Step)){Step=3}
  d[1:Step]<- cumsum(vec[1:Step])/c(1:Step)
  if(lv>Step){
    for(j in Step+1:lv){
      d[j]<- mean(vec[(j-Step+1):j])
    }
  }else{
    stop("The length of the vector cannot be longer than Step")
  }
  val<- sort(d)
  return(val)
}

#=====================================================================================================#
#Function tkvec takes quantile of tau and maps from vec to y. If vl is in between two elements in vec,
#it calculates the proportion of distance and projects into the corresponding ones in y
#' @export
tkvec<- function(vec,y,tau){
  vec<- sort(vec)
  vl<- quantile(vec,tau); vl<-unname(vl)
  id<- min(which(vec>=vl))
  rat<- abs(vec[id]-vl)/abs(vec[id]-vec[id-1]); rat=abs(rat)
  if(abs(rat) <=1){
    v = y[id]*rat + (1-rat)*y[id-1]
  }else{
    v=vec[id]
  }
  return(v)
}

#=====================================================================================================#
#' Get distribution of VaR
#'
#' This function enables one to generate the distribution of the VaR
#'
#' @param full_out full matrix of output. The rows correspond to threshold
#' values y_0 and the columns correspond to the number of MCMC draws
#' @param y_0 vector of threshold \code{y} values. It must have length equal to 
#' the number of rows of \code{full_out}
#' \code{"mode"}, or \code{"quantile"}. The default is \code{NULL} which returns the median
#' @param tau the probability tau=0.05 is the probability of VaR
#' @return val this is the distribution of the VaR
#'
#' @export

VaR<- function(full_out,y_0,tau=0.05){
  full_out<-as.matrix(full_out)
  val<- apply(full_out,2,tkvec,y_0=y_0,tau=tau)
  val<- sort(val)
  return(val)
}
