##' Solve for the SD to match rbeta to either lognormal or truncated normal
##' @param mean.val mean value 
##' @param mean.match 
##' @param sd.match
##' @param seedit
##' @param logn
##' @author Jason Cope
##' @export


opt.sd.prof<-function(mean.val,mean.match,sd.match,seedit,logn=T)
{
  obj.best<-optimize(r.beta.solve,c(0, 0.6), tol = 0.0001,mean.init=mean.val,mean.match=mean.match,sd.match=sd.match,seedit=seedit,logn=logn)[[1]]
  return(obj.best)
}