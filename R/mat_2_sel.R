##' This function calculates from the maturity values the 95%th width value needed for the logistic selectivity option in SS.
##' Use when equating maturity with selectivity	

mat_2_sel<-function(slope, Lmat50,matveccomp=seq(1,70,0.05))
{
mat.vec<-1/(1+exp(slope*(matveccomp-Lmat50)))
mat.vec[mat.vec>0.95]
width.95per<-(matveccomp[mat.vec>0.95][1]-Lmat50)
return(data.frame(Logistic_inputs=c(Lmat50,width.95per),row.names=c("Lmat","95%_width")))	
}
