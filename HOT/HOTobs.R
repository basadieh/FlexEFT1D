setwd('~/Working/FlexEFT1D/HOT')
depth_total    <- 250 #Total depth (m) of the station
dat  <-   read.csv2('HOT_obs.csv')
dat  <-   dat[-1,]
for (k in 4:ncol(dat)){
    x        <- dat[,k] 
    x        <- as.numeric(as.character(x))
    x[x<0]   <- NA
    dat[,k]  <- x
}

#Convert cast time to DOY:
dat$date  <- as.integer(as.character(dat$date))
dat$date[dat$date <0] <- NA  
dat$date  <- sprintf('%06d',dat$date)
dat$date  <- as.character(dat$date)
dat$date  <- as.Date(dat$date,'%m%d%y')
dat$DOY   <- as.numeric(strftime(dat$date, format = "%j")) 
dat$Depth <- dat$press
dat       <- dat[dat$Depth <= depth_total,]
dat0      <- dat  #Save original data

for (var in c('TIN','CHL','NPP')){
   dat <- dat0 
   if (var == 'TIN'){
      dat[,var] <- dat$nit
   
   } else if (var == 'CHL'){
   
      dat[,var] <- dat$chl
   
   } else if (var == 'NPP'){
   
      dat$d12[is.na(dat$d12)]   <- 0 
      dat[,var] <- dat$l12 + dat$d12
   
   } else{
     stop('Variable name incorrect!')
   }
   dat  <- dat[!is.na(dat[,var]),]
   dat  <- dat[!is.na(dat$DOY),  ]
   dat  <- dat[dat[,var] > 0,]
   dat1 <- dat[,c('DOY','Depth',var)]
   for (j in 1:ncol(dat1)){
       dat1[,j] <- as.numeric(dat1[,j])
   }
   if (var == 'TIN'){
      lon  <- -158
      lat  <- 22.75
      source('../getdata_station.R')
      cff  <- getdata_station('NO3',lon,lat)
      time <- cff$time 
      NO3  <- cff$data 
      
      #Correct the format:
      TIN2 <- data.frame(matrix(NA, nr = 12*nrow(NO3), nc = ncol(dat1)))
      colnames(TIN2) <- colnames(dat1)
      for (i in 1:12){
          for (j in 1:nrow(NO3)){
              TIN2[j+(i-1)*nrow(NO3),]$DOY   <- as.numeric(time[i]*30)
              TIN2[j+(i-1)*nrow(NO3),]$Depth <- abs(NO3[j,1])
              TIN2[j+(i-1)*nrow(NO3),]$TIN   <- NO3[j,i+1]
          }
      }
      TIN2 <- TIN2[TIN2$Depth <= 250,]
      dat1 <- rbind(dat1,TIN2)
   }
   #write out csv file:
   file <- paste('~/Working/FlexEFT1D/HOT/HOT','_',var,'.dat',sep='')
   write.table(dat1,file,row.names=F)
}
