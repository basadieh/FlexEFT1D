#Interpolation function for calculating MLD:
interp <- function(x,y,Xth){
       Nx <- length(x)  #Depth, must be from bottom to surface
       Ny <- length(y)  #Density difference or Kv
       stopifnot(Nx==Ny)
       z  <- NA
       
       for (i in (Nx-1):1){
           if(y[i] < Xth && y[i+1] >= Xth){
             cff <- (y[i+1]-Xth)/(y[i+1]-y[i]) 
             z   <- x[i+1]*(1-cff)+x[i]*cff
             break
           }
       }
       if (y[Nx] < Xth) z <- x[Nx] #On the surface
       return(z)
}

