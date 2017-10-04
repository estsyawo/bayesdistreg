#==============================================================================================#
#' Bayesian distribution regression
#'
#' \code{distreg} computes the distribution at a threshold value
#'
#' @param thresh threshold value that is used to binarise the continuous dependent variable
#' @param data0 original data set with the first column being the continuous dependent variable
#' @param MH metropolis-hastings algorithm to use; default:"IndepMH", alternative "RWMH"
#' @param ... any additional inputs to pass to the MH algorithm
#' @return fitob a vector of fitted values corresponding to the distribution at threshold thresh
#'
#' @examples
#' data0=DAT_ALL; qnts<-quantile(data0[,1],c(0.5,0.25,0.5,0.75,0.95));
#' distob<- distreg(qnts[1],data0)

#'
#' @export
distreg<- function(thresh,data0,MH="IndepMH",...){
  y = indicat(data0[,1],thresh) #create binary y
  data=data.frame(y,data0[,-1])#create new data set
  if(identical(MH,"IndepMH")){
    MHob=IndepMH(data = data,...)
  }else{
    MHob=RWMH(data = data,...)
  }
  fitob<-fitdist(MHob$Matpram,data)
  return(fitob)
}

#==============================================================================================#
#' Parallel compute bayesian distribution regression
#'
#' This function uses parallel computation to compute bayesian distribution regression for a given
#' vector of threshold values and a data (with first column being the continuous dependent variable)
#'
#' @param thresh vector of threshold values. ensure the min and the max generate enough zeros and ones.
#' @param data0 the original data set with a continous dependent variable in the first column
#' @param fn bayesian distribution regression function. the default is distreg provided in the package
#' @param ... any additional input parameters to pass to fn
#' @return mat a G x M matrix of output (G is the length of thresh, M is the number of draws)
#'
#' @examples
#' data0=DAT_ALL; qnts<-quantile(data0[,1],c(0.5,0.25,0.5,0.75,0.95))
#' out<- par_distreg(qnts,data0)
#'
#' @export
par_distreg<-function(thresh,data0,fn=distreg,...){ #takes a vector of threshold values
no_cores<-parallel::detectCores() - 1
c1<-parallel::makeCluster(no_cores, type = "FORK")
mat<- parallel::parSapply(c1,thresh,fn,data0=data0,...)
parallel::stopCluster(c1)
mat<- t(mat); mat<- t(apply(mat,1,sort))
return(mat)
}

#==============================================================================================#
#' Parallel compute bayesian distribution regression
#'
#' This version of the function uses the foreach package
#'
#' @param thresh vector of threshold values. ensure the min and the max generate enough zeros and ones.
#' @param data0 the original data set with a continous dependent variable in the first column
#' @param fn bayesian distribution regression function. the default is distreg provided in the package
#' @param ... any additional input parameters to pass to fn
#' @return mat a G x M matrix of output (G is the length of thresh, M is the number of draws)
#'
#' @examples
#' data0=DAT_ALL; qnts<-quantile(data0[,1],c(0.5,0.25,0.5,0.75,0.95))
#' out2<- par_distreg2(qnts,data0)
#'
#' @export
par_distreg2<-function(thresh,data0,fn=distreg,...){ #takes a vector of threshold values
  no_cores<-parallel::detectCores() - 1
  c1<-parallel::makeCluster(no_cores)
  doParallel::registerDoParallel(c1); library("foreach")
  mat=foreach(thresh = thresh, .combine = rbind,
          .export = c("thresh","data0"),.packages = "bayesdistreg")  %dopar%
    fn(thresh,data0,...)
  #mat<- parallel::parSapply(c1,thresh,fn,data0=data0,...)
  parallel::stopCluster(c1)
  mat<- t(mat); mat<- (apply(mat,1,sort))
  return(mat)
}


#==============================================================================================#
#' Parallel compute bayesian distribution regression
#'
#' \code{par_distreg3} function uses parallel computation to compute bayesian distribution regression for a given
#' vector of threshold values and a data (with first column being the continuous dependent variable). It
#' is suitable for all systems (including Windows)
#'
#' @param thresh vector of threshold values. ensure the min and the max generate enough zeros and ones.
#' @param data0 the original data set with a continous dependent variable in the first column
#' @param fn bayesian distribution regression function. the default is distreg provided in the package
#' @param ... any additional input parameters to pass to fn
#' @return mat a G x M matrix of output (G is the length of thresh, M is the number of draws)
#'
#' @examples
#' data0=DAT_ALL; qnts<-quantile(data0[,1],c(0.5,0.25,0.5,0.75,0.95))
#' out<- par_distreg3(qnts,data0)
#'
#' @export
par_distreg3<-function(thresh,data0,fn=distreg,...){ #takes a vector of threshold values
  no_cores<-parallel::detectCores() - 1
  c1<-parallel::makeCluster(no_cores, type = "PSOCK")
  mat<- parallel::parSapply(c1,thresh,fn,data0=data0,...)
  parallel::stopCluster(c1)
  mat<- t(mat); mat<- t(apply(mat,1,sort))
  return(mat)
}
