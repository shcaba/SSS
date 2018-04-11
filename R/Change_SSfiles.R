##' Change SS files
##' @param main.dir
##' @param sppname.in species name
##' @param file.make
##' @param starter.par switch to use the par file
##' @author Jason Cope
##' @export


Change.SSfiles<-function(main.dir,sppname.in,file.make=c(1,1,1,1,1),starter.par="N")

{
  for(i in 1:length(sppname.in))
  {
    if(file.make[1]==1)
    {
      if(starter.par=="N"){file.copy(from=paste(main.dir,"/baseSSfiles/starterfiles/",sppname.in[i],"_starter.ss",sep=""),to=paste(main.dir,"/",sppname.in[i],"/starter.SS",sep=""),overwrite=TRUE)}
      if(starter.par=="Y"){file.copy(from=paste(main.dir,"/baseSSfiles/starterfilespar/",sppname.in[i],"_starter.ss",sep=""),to=paste(main.dir,"/",sppname.in[i],"/starter.SS",sep=""),overwrite=TRUE)}
    }
    if(file.make[2]==1){file.copy(from=paste(main.dir,"/baseSSfiles/forecastfiles/",sppname.in[i],"_forecast.ss",sep=""),to=paste(main.dir,"/",sppname.in[i],"/forecast.SS",sep=""),overwrite=TRUE)}
    if(file.make[3]==1){file.copy(from=paste(main.dir,"/baseSSfiles/datfiles/",sppname.in[i],".dat",sep=""),to=paste(main.dir,"/",sppname.in[i],"/",sppname.in[i],".dat",sep=""),overwrite=TRUE)}
    if(file.make[4]==1){file.copy(from=paste(main.dir,"/baseSSfiles/ctlfiles/",sppname.in[i],".ctl",sep=""),to=paste(main.dir,"/",sppname.in[i],"/",sppname.in[i],".ctl",sep=""),overwrite=TRUE)}
    if(file.make[5]==1){file.copy(from=paste(main.dir,"/baseSSfiles/ss.exe",sep=""),to=paste(main.dir,"/",sppname.in[i],"/ss.exe",sep=""),overwrite=TRUE)}
  }
}