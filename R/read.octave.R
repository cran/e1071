read.octave <- function (file, quiet=FALSE) {

  nr <- 0
  nc <- 0

  if(!quiet)
    cat("Header: ")
  
  head <- scan(file=file,what=character(),nlines=4, sep=":", quiet=quiet)
  if(length(head) != 8){
    stop("Header seem to be corrupt")
  }
  for(k in 1:4){
    if(head[2*k-1] == "# rows"){
      nr <- as.integer(head[2*k])
    }
      else if(head[2*k-1] == "# columns"){
	nc <- as.integer(head[2*k])
      }
  }

  if(!quiet)
    cat("Data  : ")

  z <- scan(file=file,skip=4,quiet=quiet)
  if(length(z) != nc*nr){
    stop("Wrong number of data elements")
  }

  if((nr>1) && (nc>1)){
    if(!quiet)
      cat(paste("Matrix:", nr, "rows,", nc, "columns\n"))
    
    z<-matrix(z, nrow=nr, ncol=nc, byrow=TRUE)
  }
    else if(!quiet){
      cat("Vector:", nr*nc, "elements\n")
    }
  z
}
	      
