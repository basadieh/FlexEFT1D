require(ncdf4)
require(class)
library(fields)
# temperature data
WOATFile <- '/Users/apple/ROMS/Data/Roms_tools/WOA2009/temperature_monthly_1deg.nc'
# Nitrate data
WOANFile <- '/Users/apple/ROMS/Data/Roms_tools/WOA2009/nitrate_monthly_1deg.nc'
# Par data
parFile <- '/Users/apple/ROMS/Data/Roms_tools/SeaWifs/par_month.nc'
# Chl data
ChlFile <- '/Users/apple/ROMS/Data/Roms_tools/SeaWifs/chla_month.nc'
# Chl:C data
theta_file = '/Users/apple/Working/Global_PP/Global_size_theta.nc';
# MLD data
MLDFile <- '~/Working/Global_PP/MLD_month.nc'
# ROMS model data (one file for each month)
ROMSFile <- '/Volumes/Seagate/NPacific/NPZDSimple/npacific_avg_Y23_M'
ROMSFile <- paste(ROMSFile,1:12,'.nc',sep='')
# KD490 data
KD490File = '~/Working/Global_PP/KD490_month.nc'
POCFile   = '~/Working/Global_PP/POC_month.nc'
bathyFile <- '/Users/apple/ROMS/Data/ETOPO1_Ice_g_gmt4.grd'

#w data:
SODAFile <- '~/ROMS/Data/SODA/SODA_month.nc'
# This function does a 2-dimensional interpolation from source file to
# the desired latitude and longitude
# Improved interpolation function: directly working on the whole dataframe
Intp <- function(variable,DFin,method = 'bilinear'){
  if (class(DFin) != 'data.frame'){
    warning(paste(DFin,'must be a data.frame!'))
    DFin <- as.data.frame(DFin)
  }
  stopifnot (('Lon' %in% names(DFin)),('Lat' %in% names(DFin)))
  if (variable != 'bathymetry'){
  stopifnot ('month' %in% names(DFin))
  }
# Standardize the lon of DFin
 DFin$Lon[DFin$Lon<0] <- DFin$Lon[DFin$Lon<0]+360
 
  if (variable == 'sst'){
    sourcefile <- WOATFile
    nvar       <- 4
  } else if(variable == 'no3'){
    sourcefile <- WOANFile
    nvar       <- 4
  } else if(variable == 'par'){
    sourcefile <- parFile
    Vname      <- 'par'
    nvar       <- 1
  } else if(variable == 'chl'){
    sourcefile <- ChlFile
    nvar       <- 1
  } else if(variable == 'Chl_C'){
    sourcefile = theta_file
    nvar = 1
  } else if(variable == 'MLD'){
    sourcefile <- MLDFile
    nvar       <- 1
  } else if(variable == 'KD490'){
    sourcefile <- KD490File
    nvar       <- 1
  } else if(variable == 'POC'){
    sourcefile <- POCFile
    nvar       <- 1
  } else if(variable == 'bathymetry'){
    sourcefile <- bathyFile
    nvar       <- 1
  } else if(variable == 'w'){
    sourcefile <- SODAFile
    nvar       <- 12 
  } else if(variable == 'ssh'){
    sourcefile <- SODAFile
    nvar       <- 7 
  } else{
    stop('Variable name does not match, 
         must be one of "sst, no3, par, MLD, KD490"!')
  }
  if(file.exists(sourcefile)){
    nc       <- nc_open(sourcefile)
  }else{
    stop(paste(sourcefile,'does not exist!'))
  }
# Extract data from the nc file:
  nc  <- nc_open(sourcefile)
  v4  <- nc$var[[nvar]] # The index of data to be extracted
  lon <- v4$dim[[1]]$vals
  lon[lon<0] <- lon[lon<0]+360 #to make longitude internally consistent with different source files
  lat <- v4$dim[[2]]$vals
  var0<- ncvar_get(nc, v4)
  var0[var0==v4$missval]=NA
  nc_close(nc) 
  DFin[,variable] <- NA
  if (variable == 'bathymetry'){
    var0[var0 > 0] = 0
     d <- list(x = lon, y = lat, z = var0)  
    tmp <- DFin
# Extract Lat and Lon data from the input dataframe:
    latlon <- data.frame(lon=tmp$Lon,lat=tmp$Lat)
# Do interpolation:  
    DFin[,variable] <- interp.surface(d, latlon) 
    return(DFin)
} else{
#Loop through 12 months
  for (i in 1:12){
    tmp <- subset(DFin,month==i)
    # Extract Lat and Lon data from the input dataframe:
    latlon <- data.frame(lon=tmp$Lon,lat=tmp$Lat)

    if(nrow(tmp) < 1) {
      next
    } else{
    if (method == 'bilinear'){
# Check the dimensions of target data
    if (v4$ndims == 3){
      d <- list(x = lon, y = lat, z = var0[,,i])
    }else if(v4$ndims ==4){
      # Extract the surface values
      d <- list(x = lon, y = lat, z = var0[,,1,i])
    }else{
      stop('Data dimensions wrong! Please check!')
    }
# Using the bilinear interpolation method
# This method is more reliable than the nearest neibour method 
# given that the value to be estimated has to be within a grid 
# where four corners must have values


# Do interpolation:  

    DFin[DFin$month==i,][,variable] <- interp.surface(d, latlon) 

    } else if(method == 'nearest'){

    pts  <- expand.grid(lon=lon, lat=lat)
    # Check the dimensions of target data
    if (v4$ndims == 3){
      d <- var0[,,i]
    }else if(v4$ndims ==4){
      # Extract the surface values
      d <- var0[,,1,i]
    }else{
      stop('Data dimensions wrong! Please check!')
    }

    DFin[DFin$month==i,][,variable] <- 
       as.numeric(as.character(knn(pts[!is.na(d),], 
                            latlon, d[!is.na(d)])))
   }
  }
 }
#End loop
  return(DFin)
}
}

# A useful function reads from a nc file and generates an list that contains lon,lat and the real data (3D array)
readnc <- function(variable){
  if (variable == 'sst'){
    sourcefile <- WOATFile
    Vname      <- 't_an'
    nvar       <- 4
  } else if(variable == 'Aks'){
    sourcefile <- ROMSFile[1]
    Vname      <- 'AKs'
    lonname    <- 'lon_rho'
    latname    <- 'lat_rho'
  } else if(variable == 'Chl_C'){
    sourcefile <- theta_file
    Vname      <- 'Chl_C'
    nvar       <- 1
  } else if(variable == 'w_SODA'){
    sourcefile <- SODAFile
    nvar       <- 12 
  } else if(variable == 'w_ROMS'){
    sourcefile <- ROMSFile[1]
    Vname      <- 'w'
    lonname    <- 'lon_rho'
    latname    <- 'lat_rho'

  } else if(variable == 'ssh'){
    sourcefile <- SODAFile
    nvar       <- 7 
  } else if(variable == 'no3'){
    sourcefile <- WOANFile
    nvar       <- 4
  } else if(variable == 'par'){
    sourcefile <- parFile
    Vname      <- variable
    nvar       <- 1
  } else if(variable == 'chl'){
    sourcefile <- ChlFile
    nvar       <- 1
  } else if(variable == 'MLD'){
    sourcefile <- MLDFile
    nvar       <- 1
  } else if(variable == 'KD490'){
    sourcefile <- KD490File
    nvar       <- 1
  } else if(variable == 'POC'){
    sourcefile <- POCFile
    nvar       <- 1
  } else{
    stop('Variable name does not match, 
         must be one of "sst, no3, par, MLD, KD490"!')
  }
  if(file.exists(sourcefile)){
    nc       <- nc_open(sourcefile)
  }else{
    stop(paste(sourcefile,'does not exist!'))
  }
  if (!exists('nvar')) {
      nvar <- which(names(nc$var) == Vname)
  }
  v4  <- nc$var[[nvar]] # The index of data to be extracted

  if (!(variable %in% c('Aks','w_ROMS'))){
      lon <- v4$dim[[1]]$vals
      lat <- v4$dim[[2]]$vals

  } else{
      lon <- ncvar_get(nc,lonname)
      lat <- ncvar_get(nc,latname)
      lon <- lon[,1]
      lat <- lat[1,]
  }
  lon[lon<0]             <- lon[lon<0]+360 #to make longitude internally consistent with different source files
  var0                   <- ncvar_get(nc, v4)
  var0[var0==v4$missval] <- NA

  if (length(dim(var0)) == 3) {

    time  <- v4$dim[[3]]$vals
    return(list(lon=lon,lat=lat,time=time,data=var0))

  } else if (length(dim(var0)) == 4) {

    depth <- v4$dim[[3]]$vals

    if (sourcefile == SODAFile){
       depth <- ncvar_get(nc,'Depth')
    }
    time  <- v4$dim[[4]]$vals

    # Need to convert depth for ROMS output
    if( (variable == 'Aks') | (variable == 'w_ROMS') ) {
       hc = ncvar_get(nc,'hc') 
     Cs_r = ncvar_get(nc,'Cs_r')
     sc_r = ncvar_get(nc,'sc_r')
       h  = ncvar_get(nc,'h') 
# Convert ROMS coordinate to real depth, assuming zeta = 0:
    depth = array(NA,c( dim(h), dim(Cs_r) ))
    for (k in 1:length(Cs_r) ){
        depth[,,k] = (sc_r[k] - Cs_r[k])*hc + Cs_r[k]*h
        }
   }
    
    if (variable %in% c('Aks','w_ROMS')){
       nc_close(nc) 
    
       # Estimate time from multiple nc files
       time1         <- numeric(length(time) * 12)
       time1[1:30]   <- time
       time          <- time1
       
       # Get data from multiple nc files
       var1          <- array(NA,dim=c(dim(var0)[1:3],dim(var0)[4]*12))
       var1[,,,1:30] <- var0
       var0          <- var1
       rm(var1)
       for (mo in 2:12){
           
            nc <- nc_open(ROMSFile[mo])
           v4  <- nc$var[[nvar]] # The index of data to be extracted
          var1 <- ncvar_get(nc, v4)
          var1[var1==v4$missval] <- NA
            time[ ((mo-1)*30+1):(mo*30)] <- v4$dim[[4]]$vals
          var0[,,,((mo-1)*30+1):(mo*30)] <- var1
          nc_close(nc)
       }
    }
    return(list(lon=lon,lat=lat,time=time,depth=depth,data=var0))
  } else {
    stop('Data dimension incorrect!')
  }
  if (!(variable %in% c('Aks','w_ROMS'))) {nc_close(nc) }
}



# Read lon and lat from grid file:
read_grid <- function(nameit){
  grdname = paste('~/Roms_tools/Run/ROMS_FILES/Output/',
                  nameit,'-grid.nc',sep='')
  if(file.exists(grdname)){
    ng=nc_open(grdname)
  }else{
    stop(paste(grdname,'does not exist!'))
  }
  lon=ng$var$lon_rho
  lat=ng$var$lat_rho
  lon=ncvar_get(ng,lon)
  lon=lon[,1]
#  lon[lon<0] <- lon[lon<0]+360
  lat=ncvar_get(ng,lat)
  lat=lat[1,]
  nc_close(ng)
  return(list(lon=lon,lat=lat))
}

#the following function fill_grid generates a 3D array to fill in the data of a regular rectangle surface
fill_grid <- function(grid,varname,method='bilinear'){
  #lon and lat of MLD (Training data)
  grid$lon[grid$lon > 0] <- grid$lon[grid$lon > 0] - 360
  MLD = readnc(varname)
  LON = MLD$lon; LAT=MLD$lat
  LON[LON>0]=LON[LON>0]-360
  dat = array(NA,c(length(grid$lon),length(grid$lat),12)) 
  if (method == 'nearest'){
    # Test data:
    latlon = expand.grid(lon=grid$lon,lat=grid$lat)
    pts = expand.grid(lon=LON, lat=LAT)
    for (i in 1:12){
        var = MLD$data[,,i]
        #tmp is a vector
        tmp = as.numeric(as.character(knn(pts[!is.na(var),],latlon, 
                                          var[!is.na(var)])))
        #Convert tmp to a matrix
        dat[,,i] = matrix(tmp,nr=length(grid$lon),nc=length(grid$lat))
    }
  } else if(method == 'bilinear'){
# This method much faster
    Z = MLD$data
    for (i in 1:12){
      if(length(dim(Z)) == 3) {
       d=list(x=LON,y=LAT,z=Z[,,i])
      }else if (length(dim(Z)) == 4){
       d=list(x=LON,y=LAT,z=Z[,,1,i])    
      }
      gridlist=list(x=grid$lon,y=grid$lat)
      t <- interp.surface.grid(d,gridlist)
      dat[,,i]=t$z
    }
  } else{
    stop('Selected method is incorrect!')
  }
    return(dat)
}

