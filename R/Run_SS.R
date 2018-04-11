##' Run SS function
##' @param path location of the executable
##' @param ss.exe executable name for ss (e.g., "ss")
##' @param ss.cmd command line input used when running SS
##' @author Jason Cope
##' @export

RUN.SS<-function(path, ss.exe="ss",ss.cmd=" -nohess -nox")
{
  navigate <- paste("cd ",path,sep="")
  command <- paste(navigate," & ",ss.exe,ss.cmd,sep="")
  shell(command,invisible=TRUE,translate=TRUE)
}