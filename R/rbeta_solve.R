##' Objective function used to match SD from beta to either lognormal or truncated normal
##' @param sd.init
##' @param mean.init
##' @param mean.match
##' @param sd.match
##' @param seedit
##' @param logn
##' @author Jason Cope
##' @export

#Objective function used to match SD from beta to either lognormal or truncated normal
r.beta.solve<-function(sd.init,mean.init,mean.match,sd.match,seedit,logn=T)
{
  set.seed(seedit)
  tbeta<-rbeta.ab(100000,mean.init,sd.init,0.99,0.01)
  if(logn==T){log_ns<-rlnorm(1000000,log(mean.match),sd.match)}
  if(logn==F){log_ns<-rtnorm(1000000,mean.match,sd.match)}
  mean.diff<-abs(mean(log_ns)-mean(tbeta))
  sd.diff<-abs(sd(log_ns)-sd(tbeta))
  med.diff<-abs(median(log_ns)-median(tbeta))
  obj<--(mean.diff+sd.diff+med.diff)
  return(obj)
}