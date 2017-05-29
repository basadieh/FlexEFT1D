setwd('~/Working/FlexEFT1D')

#Primary production data:
PP       <- read.csv('Obs_data/PP.csv')
PP$Date  <- as.Date(PP$Date,'%m/%d/%Y')
PP$DOY   <- as.numeric(strftime(PP$Date, format = "%j")) 
PP$Depth <- PP[,'Inc..depth']
PP$PP1   <- PP$PP.1
PP$PP2   <- as.numeric(as.character(PP$PP.2))
PP$PP    <- sapply(1:nrow(PP),function(i)mean(c(PP[i,]$PP1,PP[i,]$PP2),na.rm=T))
S1PP     <- PP[PP$Stn.=='S1',c('DOY','Depth','PP')]
K2PP     <- PP[PP$Stn.=='K2',c('DOY','Depth','PP')]
write.table(S1PP,'S1/S1_PP.dat',row.names=F)
write.table(K2PP,'K2/K2_PP.dat',row.names=F)

#Size-fraction Chl data:
pigment      <- read.csv('Obs_data/Pigment.csv')
pigment$DATE <- as.Date(pigment$DATE,'%m/%d/%Y')
pigment$DOY  <- as.numeric(strftime(pigment$DATE, format='%j'))
pigment$Depth<- pigment$CTDDPT
pigment      <- pigment[(pigment$Depth >=0) & (!is.na(pigment$Depth)), ]

for (i in c('SIZE10','SIZE3','SIZE1','SIZE_L1')){
  pigment[,i][pigment[,i] < 0] <- NA
}

pig2   <- pigment[(!is.na(pigment$SIZE10)) | (!is.na(pigment$SIZE3)) 
                  |(!is.na(pigment$SIZE1)) | (!is.na(pigment$SIZE_L1)),]

S1size <- pig2[pig2$STNNBR=='S01',
                  c('DOY','Depth','SIZE10','SIZE3','SIZE1','SIZE_L1')]


K2size <- pig2[pig2$STNNBR=='K02',
               c('DOY','Depth','SIZE10','SIZE3','SIZE1','SIZE_L1')]

#write.csv(S1size,'S1/S1_size.csv',row.names=F)
#for (i in 3:ncol(S1size)){
#  z          <- S1size[,i]
#  S1size[,i] <- sapply(1:length(z),function(x)max(z[x],1e-3))
#}
#for (i in 3:ncol(K2size)){
#  z          <- K2size[,i]
#  K2size[,i] <- sapply(1:length(z),function(x)max(z[x],1e-3))
#}

write.table(K2size,'K2/K2_size.dat',row.names=F)
write.table(S1size,'S1/S1_size.dat',row.names=F)

####==================================================================================
#Chl, NO3 data for K2,S1
get_CHLNO3 <- function(stn){  
  datfile <- paste('~/working/FlexEFT1D/Obs_data/',stn,'_obs.csv',sep='')
  
  if ( file.exists(datfile) ){
        dat = read.csv(datfile)
      if (ncol(dat) <= 1){
        dat = read.csv2(datfile)
      }
  } else{
      stop(paste('Observation file',datfile,' unavailable!',sep=''))
  }
  dat$Lon   = as.numeric(as.character(dat$Longitude))
  dat$Lat   = as.numeric(as.character(dat$Latitude))
  dat$Date  = as.Date(dat$Date,'%m/%d/%Y')
  dat$DOY   = as.numeric(strftime(dat$Date, format = "%j"))    #Date of the year
  dat$Depth = as.numeric(as.character(dat$Depth)) 
  #Only need data above 250 m
  dat = dat[dat$Depth <= 250,]
  
  #Clean up nitrate data:
  dat$Nitrate=as.numeric(as.character(dat$Nitrate))
  dat$Nitrate[dat$Nitrate<0]=NA
  dat$Nitrate1=as.numeric(as.character(dat$Nitrate1))
  dat$Nitrate1[dat$Nitrate1<0]=NA
  
  dat <- dat[!is.na(dat$Nitrate),]
  #Clean up nitrite data:
  dat$Nitrite=as.numeric(as.character(dat$Nitrite))
  dat$Nitrite[dat$Nitrite<0]=NA
  
  dat$Nitrite1=as.numeric(as.character(dat$Nitrite1))
  dat$Nitrite1[dat$Nitrite1<0]=NA
  
  #Clean up NH4 data:
  dat$NH4=as.numeric(as.character(dat$NH4))
  dat$NH4[dat$NH4<0]=NA
  dat$NH41=as.numeric(as.character(dat$NH41))
  dat$NH41[dat$NH41<0]=NA
  
  
  dat$NO3 <- sapply(1:nrow(dat),function(i)mean(c(dat$Nitrate[i], dat$Nitrate1[i]),na.rm=T)) 
  dat$NO2 <- sapply(1:nrow(dat),function(i)mean(c(dat$Nitrite[i], dat$Nitrite1[i]),na.rm=T)) 
  dat$Nh4 <- sapply(1:nrow(dat),function(i)mean(c(dat$NH4[i], dat$NH41[i]),na.rm=T)) 
  
  #Add up
  dat$TIN <- sapply(1:nrow(dat),
              function(i)sum(c(dat$NO3[i], dat$NO2[i], dat$Nh4[i]),na.rm=T)) 
  
  #Clean up TIN data
  dat = dat[dat$TIN > 0,]
  TIN = dat[,c('DOY','Depth','TIN')]
  
  #Extract WOA13 data:
  
  depth_total    <- 250 #Total depth (m) of the station
  
  if (stn == 'S1'){
      lon  = 145
      lat  = 30
  } else if(stn == 'K2'){
      lon  = 160 
      lat  = 47
  } else if(stn == 'HOT'){
      lon  = -158
      lat  = 22.75
  }
  
  source('~/Working/FlexEFT1D/Rscripts/getdata_station.R')
  cff  <- getdata_station('NO3',lon,lat)
  time <- cff$time 
  NO3  <- cff$data 
  
  #Correct the format:
  TIN2 <- data.frame(matrix(NA, nr = 12*nrow(NO3), nc = ncol(TIN)))
  colnames(TIN2) <- colnames(TIN)
  for (i in 1:12){
      for (j in 1:nrow(NO3)){
          TIN2[j+(i-1)*nrow(NO3),]$DOY   <- time[i]*30
          TIN2[j+(i-1)*nrow(NO3),]$Depth <- abs(NO3[j,1])
          TIN2[j+(i-1)*nrow(NO3),]$TIN   <- NO3[j,i+1]
      }
  }
  TIN <- rbind(TIN,TIN2)
  for (i in 1:ncol(TIN)){
      TIN[,i] <- as.numeric(as.character(TIN[,i]))
  }
  
  TIN <- TIN[TIN$Depth <= depth_total, ]
  #write out dat file:
  filedir <- paste('~/Working/FlexEFT1D/',stn,'/',stn,'_TIN.dat',sep='')
  write.table(TIN,file = filedir,row.names=F)
  
  #-------------------------
  #Set up CHL data
  if ( file.exists(datfile) ){
        dat <- read.csv(datfile)
      if (ncol(dat) <= 1){
        dat <- read.csv2(datfile)
      }
  } else{
        stop(paste('Observation file',datfile,' unavailable!',sep=''))
  }
  
  dat$Lon   = as.numeric(as.character(dat$Longitude))
  dat$Lat   = as.numeric(as.character(dat$Latitude))
  dat$Date  = as.Date(dat$Date,'%m/%d/%Y')
  dat$DOY   = as.numeric(strftime(dat$Date, format = "%j"))    #Date of the year
  dat$Depth = as.numeric(as.character(dat$Depth)) 
  #Only need data above 250 m
  dat       = dat[dat$Depth <= depth_total,]
  
  dat$CHL   = as.numeric(as.character(dat$CHLWEL))
  dat$CHL[dat$CHL<0]=NA
  dat$CHL1  = as.numeric(as.character(dat$CHLWEL1))
  dat$CHL1[dat$CHL1<0]=NA
  dat$Chl   = sapply(1:nrow(dat),function(i)mean(c(dat$CHL[i], dat$CHL1[i]),na.rm=T)) 
  dat       = dat[dat$Chl >= 0,]
  CHL       = dat[,c('DOY','Depth','Chl')]
  CHL       = CHL[!is.na(CHL$Chl),]
  
  # Extract data from MODIS:
  #
  #Open the file
  file <- paste('~/Working/FlexEFT1D/Obs_data/',stn,'_Chl_modis.csv',sep='')
  dat  <- read.csv(file)
  if (ncol(dat) <= 1){
     dat <- read.csv2(file)
  }
  #Calculate DOY:
  dat$Date  <- as.Date(dat$Date,'%Y/%m/%d')
  #Date of the year
  dat$DOY   <- as.numeric(strftime(dat$Date, format = "%j"))  
  dat$Depth <- 3
  #Remove NA data:
  dat$Chl   <- dat$MeanChl
  dat       <- dat[!is.na(dat$Chl),]
  dat       <- dat[,c('DOY','Depth','Chl')]
  dat       <- dat[dat$DOY <= 360,]
  #write out csv file:
  CHL       <- rbind(CHL,dat)
  
  filedir <- paste('~/Working/FlexEFT1D/',stn,'/',stn,'_CHL.dat',sep='')
  write.table(CHL,filedir,row.names=F)
}

#Test likelihood function:
p <- function(n,sigma,SS){
        return(-n*log(sigma*sqrt(2*pi)) - SS/2/sigma**2)
     }
