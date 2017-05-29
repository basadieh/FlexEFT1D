#Write a 1D code using R without physical computation, read in external parameters
#and compute biological parameters only
setwd('~/Working/FlexEFT1D/S1')

#Compile external fortran file:
#system('ifort -c FlexEFT.f90')
#system('ifort -shared -o FlexEFT.so FlexEFT.o')

#1. Set up model run duration, Lon and Lat, model time step, total depth, and number of depths
Stn <- 'S1'
lon <- 145
lat <- 30
  
source('~/Working/Global_PP/Interpolate_WOA.R')
#For each depth, get the profile at the targeted coordinates
Temp_data <- readnc('sst')
 Par_data <- readnc('par')
 NO3_data <- readnc('no3')
 Aks_data <- readnc('Aks')
#Average Aks to monthly intervals:
m_per_s   <- 1/(3600*24*30)
time      <- Aks_data$time*m_per_s
time      <- time%%12
month     <- numeric(12)
Aks       <- Aks_data$data
Aks       <- array(NA,c(dim(Aks)[1],dim(Aks)[2],dim(Aks)[3],12))
for (i in 1:12){
 month[i]  <- mean(time[c((3*(i-1)+1):(3*i))])
 data      <- Aks_data$data[,,,c((3*(i-1)+1):(3*i))]
 Aks[,,,i] <- apply(data,1:3,function(x)mean(x,na.rm=T))
}
Aks[Aks <= 1e-4] <- 0.0001
Aks_data$time  <- month
Aks_data$data  <- Aks
Temp_data$time <- month
NO3_data$time  <- month
save(Temp_data,Par_data,Aks_data,NO3_data,file = 'Input.Rdata')

load('Input.Rdata')
getdata_station <- function(varname,lon,lat){
    #First, get the total dataset
    if (varname == 'temp'){
      data_s <- Temp_data
    } else if (varname == 'par'){
      data_s <- Par_data
    } else if (varname == 'Aks'){
      data_s <- Aks_data
    } else if (varname == 'NO3'){
      data_s <- NO3_data
    } else{
      stop('Variable name incorrect!')
    }
    Lon_s  <- data_s$lon
    Lat_s  <- data_s$lat
    Lon_s[Lon_s>0] <- Lon_s[Lon_s>0]-360
    # Obtain the data:
    Z      <- data_s$data

    # Determine the time dimension of the data
    time       <- data_s$time
    lon[lon>0] <- lon[lon>0]-360
    gridlist   <- list(x=lon,y=lat)

    if(length(dim(Z)) == 3) {
       dat  <- data.frame(matrix(NA, nr=1, nc=length(time)+1 ) )
      names(dat) <- c('Depth', paste('M',1:12,sep=''))
       dat[,1]   <- 0

      for (i in 1:length(time) ) {
          d    <- list(x=Lon_s, y=Lat_s, z=Z[,,i])
        dat[1,i+1] <- interp.surface.grid(d,gridlist)$z
      }

     } else if (length(dim(Z)) == 4){
      depth <- data_s$depth
      if (varname == 'Aks'){
        N = dim(depth)[3]
       dat= numeric(N)
        for (k in 1:N){
            d      <- list(x=Lon_s, y=Lat_s, z=depth[,,k] )    
            dat[k] <- interp.surface.grid(d,gridlist)$z 
        }
       depth = abs(dat)
      } 
       dat  <- data.frame(matrix(NA, nr=length(depth), nc=length(time)+1 ) )

      names(dat) <- c('Depth', paste('M',1:length(time),sep=''))
       dat[,1]   <- -depth
      for (i in 1:length(depth)) {
          for (j in 1:length(time) ){
              d <- list(x=Lon_s, y=Lat_s, z=Z[,,i,j])    
            dat[i,j+1] <- interp.surface.grid(d,gridlist)$z
          }
      }
     } else{
     stop('Data dimension incorrect!')
     }
     return(dat)
}


#Write into data files:
for (var in c('temp','par','Aks','NO3')){
    data     <-  getdata_station(var,lon,lat) 
    if (var == 'NO3'){data = data[,1:2]}
    if ((var != 'Aks') && (var != 'par')){
       data = data[nrow(data):1,]
    }
    outfile  <-  paste(Stn,'_',var,'.dat',sep='')
    write.table(data,file=outfile,append = F,row.names=FALSE,col.names=TRUE) 
}

#Plot these variables:

#Plot NO3:
var      <- 'NO3'
data     <-  getdata_station(var,lon,lat) 
library(plot3D)
depth <- data$Depth
month <- seq(0.5,11.5,1)
data1 <- as.matrix(t(data[nrow(data):1,2:ncol(data)]))
image2D(data1, x=month, y=rev(depth),resfac=5,yaxt='n',
            xlab="Month", ylab="Depth (m)")
axis(2, at = seq(0,min(depth),by=-50)) 

plot(month,data1[,ncol(data1)],type='b')


get_external_profile <- function(varname,lon,lat,doy){
    #First, get the total dataset
    if (varname == 'temp'){
      data_s <- Temp_data
    } else if (varname == 'par'){
      data_s <- Par_data
    } else if (varname == 'Aks'){
      data_s <- Aks_data
    } else{
      stop('Variable name incorrect!')
    }
    Lon_s  <- data_s$lon
    Lat_s  <- data_s$lat
    Lon_s[Lon_s>0] <- Lon_s[Lon_s>0]-360
    Z      <- data_s$data
    doy    <- as.integer(doy)
    month  <- (doy-1)%/%30 + 1
    day    <- doy-30*(month - 1)
    lon[lon>0] <- lon[lon>0]-360
    gridlist <- list(x=lon,y=lat)

    if (day >= 15){
         month2 <- month + 1
       if (month2 > 12) month2 <- 1
    }  else{
         month2 <- month - 1
       if (month2 < 0) month2 <- 12
     }

    if(length(dim(Z)) == 3) {

      d1 <- list(x=Lon_s, y=Lat_s, z=Z[,,month])
     t1  <- interp.surface.grid(d1,gridlist)$z
      d2 <- list(x=Lon_s, y=Lat_s, z=Z[,,month2])
     t2  <- interp.surface.grid(d2,gridlist)$z
     dat <- (t1*(31-abs(day-15)) + t2*abs(day-15))/31
     
     } else if (length(dim(Z)) == 4){
     depth <- data_s$depth
     dat   <- list(depth = depth, data = numeric(length(depth)))
      for (i in 1:length(depth)) {
          d1 <- list(x=Lon_s, y=Lat_s, z=Z[,,i,month])    
          d2 <- list(x=Lon_s, y=Lat_s, z=Z[,,i,month2])    
          t1 <- interp.surface.grid(d1,gridlist)$z
          t2 <- interp.surface.grid(d2,gridlist)$z
      dat$data[i] <- (t1*(31-abs(day-15)) + t2*abs(day-15))/31
      }
     } else{
     stop('Data dimension incorrect!')
     }
     return(dat)
}

#Test:
system.time(get_external_profile('temp',Lon,Lat,20))
system.time(get_external_profile('temp',Lon,Lat,350))
system.time(get_external_profile('Aks',Lon,Lat,350))
system.time(get_external_profile('par',Lon,Lat,75))


#Wrap up to calculate the new set of state variables:
biorhs <- function(temp,par,dtdays,NVAR,Vars){
          dyn.load('FlexEFT.so')
          newVAR <- .Fortran('FlexEFT',tC = as.numeric(temp),
                      par    = as.numeric(par),
                      dtdays = as.numeric(dtdays),NVar = as.integer(NVAR),
                      Vars   = as.numeric(Vars))$Vars
          return(newVAR)
          }
#7. If savetime, save data into the file.


