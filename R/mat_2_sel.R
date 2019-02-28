##' This function calculates from the maturity values the 95%th width value needed for the logistic selectivity option in SS.
##' Use when equating maturity with selectivity	

mat_2_sel<-function(slope, Lmat50,matveccomp=seq(1,70,0.05))
{
mat.vec<-1/(1+exp(slope*(x-intercept)))
mat.vec[mat.vec>0.95]
width.95per<-(x[mat.vec>0.95][1]-intercept)
return(data.frame(Logistic_inputs=c(Lmat50,width.95per),row.names=c("Lmat","95%_width")))	
}
