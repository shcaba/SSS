##' This function grabs a table of values out of the report file	

  matchfun2 <- function(string1,adjust1,string2,adjust2,cols="nonblank",matchcol1=1,matchcol2=1,
    objmatch=rawrep,objsubset=rawrep,substr1=TRUE,substr2=TRUE,header=FALSE)
  	{
  	  	# return a subset of values from the report file (or other file)
  	  	# subset is defined by character strings at the start and end, with integer
  	  	# adjustments of the number of lines to above/below the two strings
  	  	line1 <- match(string1,
  	  	               if(substr1){
  	  	                 substring(objmatch[,matchcol1],1,nchar(string1))
  	  	               }else{
  	  	                 objmatch[,matchcol1]
  	  	               })
  	  	line2 <- match(string2,
  	  	               if(substr2){
  	  	                 substring(objmatch[,matchcol2],1,nchar(string2))
  	  	               }else{
  	  	                 objmatch[,matchcol2]
  	  	               })
  	  	if(is.na(line1) | is.na(line2)) return("absent")
		
  	  	if(is.numeric(cols))    out <- objsubset[(line1+adjust1):(line2+adjust2),cols]
  	  	if(cols[1]=="all")      out <- objsubset[(line1+adjust1):(line2+adjust2),]
  	  	if(cols[1]=="nonblank"){
  	  	  # returns only columns that contain at least one non-empty value
  	  	  out <- objsubset[(line1+adjust1):(line2+adjust2),]
  	  	  out <- out[,apply(out,2,emptytest) < 1]
  	  	}
  	  	if(header && nrow(out)>0){
  	  	  out[1,out[1,]==""] <- "NoName"
  	  	  names(out) <- out[1,]
  	  	  out <- out[-1,]
  	  	}
  	  	return(out)
  	}