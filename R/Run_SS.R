##' Run SS function
##' @param path location of the executable
##' @param ss.exe executable name for ss (e.g., "ss")
##' @param ss.cmd command line input used when running SS
##' @author Jason Cope
##' @export

RUN.SS<-function(path,ss.exe="ss",ss.cmd=" -nohess -nox",OS.in=NA){ 
  navigate <- paste("cd ", path, sep="") 
if(OS.in=="Windows") 
  {
    command <- paste0(navigate," & ", ss.exe, ss.cmd) 
    shell(command, invisible=TRUE, translate=TRUE)
  } 
if(OS.in=="Mac")  
  {
    
    command <- c(paste("cd", path), "chmod +x ./ss_mac","./ss_mac") 
    system(paste(command, collapse=";"),invisible=TRUE)
    
    #command <- paste0(path,"/./ss_mac", ss.cmd) 
    #system(command, invisible=TRUE)
  } 
if(OS.in=="Linux") 
  {
    command <- c(paste("cd", path), "chmod +x ./ss_linux","./ss_linux") 
    system(paste(command, collapse=";"), invisible=TRUE)
  } 
  
}  



# RUN.SS<-function(path, ss.exe="ss",ss.cmd=" -nohess -nox")
# {
#   navigate <- paste("cd ",path,sep="")
#   command <- paste(navigate," & ",ss.exe,ss.cmd,sep="")
#   shell(command,invisible=TRUE,translate=TRUE) 	
# }