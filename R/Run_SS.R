##' Run SS function
##' @param path
##' @param ss.exe
##' @param ss.cmd
##' @author Jason Cope
##' @export

RUN.SS<-function(path, ss.exe="ss",ss.cmd=" -nohess -nox")
{
  navigate <- paste("cd ",path,sep="")
  command <- paste(navigate," & ",ss.exe,ss.cmd,sep="")
  shell(command,invisible=TRUE,translate=TRUE)
}