#Read in external parameters
#and compute biological parameters only
setwd('~/Working/FlexEFT1D')
library(plot3D)
library(oce)
source('Rscripts/getdata_station.R')
source('Rscripts/interp.R')
source('Rscripts/Hz.R')

WOA13NFile <- '~/ROMS/Data/WOA13/WOA13_NO3.Rdata'
load(WOA13NFile)

WOA13TFile <- '~/ROMS/Data/WOA13/WOA13_Temp.Rdata'
load(WOA13TFile)

#Compile external fortran file:
#system('ifort -c FlexEFT.f90')
#system('ifort -shared -o FlexEFT.so FlexEFT.o')

Stn_name <- 'S1'
depth    <- 250 #Total depth (m) of the station

if (Stn_name == 'S1'){
    stn_lon  = 145
    stn_lat  = 30
} else if(Stn_name == 'K2'){
    stn_lon  = 160 
    stn_lat  = 47
} else if(Stn_name == 'HOT'){
    stn_lon  = -158
    stn_lat  = 22.75
}

#For each depth, get the profile at the targeted coordinates
 Temp_data <- woa13temp
 Par_data  <- readnc('par')
 NO3_data  <- woa13no3

 Aks_data  <- readnc('Aks')       #From my ROMS output, unit: m2/s
wSODA_data <- readnc('w_SODA')    #unit: m/s
wROMS_data <- readnc('w_ROMS')    #unit: m/s

#Correct timing for Aks and wROMS
m_per_s   <- 1/(3600*24*30)
time      <- Aks_data$time*m_per_s
time      <- time%%12
  Aks_data$time  <- time
wROMS_data$time  <- time

  #Calculate MLD based on observed profiles:
MLD_obs <- function(Stn_name){
   datfile <- paste('~/working/FlexEFT1D/Obs_data/',Stn_name,
                    '_obs.csv',sep='')
   
   if ( file.exists(datfile) ){
         dat = read.csv(datfile)
       if (ncol(dat) <= 1){
         dat = read.csv2(datfile)
       }
   } else{
       stop(paste('Observation file',datfile,' unavailable!',sep=''))
   }
   dat$Date   = as.Date(dat$Date,'%m/%d/%Y')
   dat$DOY    = as.numeric(strftime(dat$Date, format = "%j"))    #Date of the year
   dat$Depth  = as.numeric(as.character(dat$Depth)) 

   #Only need data above 250 m
   dat        = dat[dat$Depth <= 250,]
   dat$Sigma  = dat$Sigma_Theta
   dat$Sigma[dat$Sigma < 0] = NA
   Date_uni   = unique(dat[,c('Station','Date','DOY')])
   N          = nrow(Date_uni)
   MLD_obs    = data.frame(matrix(NA,nr=N,nc=4))
   names(MLD_obs) = c('Station','Date','DOY','MLD')
   MLD_obs$Station= as.character(MLD_obs$Station)
   MLD_obs$Date   = as.Date(MLD_obs$Date)
   MLD_obs$DOY    = as.integer(MLD_obs$DOY)
   sigma_th   = 0.125  #Threshold of density for MLD
   for (i in 1:N){
       k      = (dat$Station==Date_uni[i,'Station']) &
                (dat$Date   ==Date_uni[i,'Date'])    &
                (dat$DOY    ==Date_uni[i,'DOY'])
       DOY    = dat[k,'DOY']
       stopifnot (sd(DOY) == 0) 
       MLD_obs$DOY[i]     = DOY[1]
       MLD_obs$Station[i] = as.character(dat[k,'Station'][1])
       MLD_obs$Date[i]    = dat[k,'Date'][1]
       dep    = dat[k,'Depth']
       sigma  = dat[k,'Sigma']
       Z      = data.frame(Depth = -abs(dep), density = sigma)
       Z      = na.omit(Z)
       Z      = Z[order(Z$Depth),]
       Z$dif  = Z$density-Z$density[nrow(Z)]
       MLD_found = FALSE
       for (j in (nrow(Z)-1):1){
           if ((Z$dif[j] >= sigma_th) && (Z$dif[j+1] < sigma_th)){
              #Find the layer satisfying the MLD threshold
              cff <- (Z$dif[j]-sigma_th)/(Z$dif[j]-Z$dif[j+1]) 
              MLD <- Z$Depth[j]*(1-cff)+Z$Depth[j+1]*cff
              MLD_found <- TRUE
           } 
           if (MLD_found)  break
       }
       MLD_obs$MLD[i] = MLD
   }
   return(MLD_obs)
}

#Write into data files:
for (var in c('temp','par','Aks','NO3','wROMS','wSODA')){
    cff  = getdata_station(var,stn_lon,stn_lat)
    time = cff$time 
    data = cff$data 
    if (!(var %in% c('Aks','wROMS', 'par')) ){
       data = data[nrow(data):1,]          #Bottom layer first
    }

    outfile  <-  paste('~/Working/FlexEFT1D/',
                        Stn_name,'/',Stn_name,'_',var,'.dat',sep='')
    write.table(data,file=outfile,  
                  append = F,row.names=FALSE,col.names=TRUE) 
    data      <-  read.table(outfile,header=T)
    #Write out timefile:
    timefile  <-  paste('~/Working/FlexEFT1D/',
                         Stn_name,'/',Stn_name,'_',var,'_time.dat',sep='')
    write.table(time,file=timefile,  
                  append = F,row.names=FALSE,col.names=TRUE)  
    time      <-  read.table(timefile,header=T)
}

plot_forcing <- function(Stn_name,var,useTaketo=FALSE){
    outfile  =  paste('~/Working/FlexEFT1D/',
                      Stn_name,'/',Stn_name,'_',var,'.dat',sep='')
    #Write out timefile:
    timefile =  paste('~/Working/FlexEFT1D/',
                      Stn_name,'/',Stn_name,'_',var,'_time.dat',sep='')

    pdffile  = paste('~/Working/FlexEFT1D/',
                     Stn_name,'/',Stn_name,'_',var,'.pdf',sep='')

    if (useTaketo && var == 'Aks'){
     timefile= paste('~/Working/FlexEFT1D/',Stn_name,'/',
                          Stn_name,'_Aks_timeT.dat',sep='')

     outfile = paste('~/Working/FlexEFT1D/',Stn_name,'/',Stn_name,'_AksT.dat',sep='')
     pdffile = paste('~/Working/FlexEFT1D/',
                     Stn_name,'/',Stn_name,'_',var,'T.pdf',sep='')

    }
    time     =  read.table(timefile,header=T)
    data     =  read.table(outfile, header=T)

       pdf(pdffile,width=5,height=5,paper='a4')
          par(font.lab  = 1,
              family    = "serif",
              mar       = c(4,4,4,4),
              mgp       = c(2.2,1,0))
          time  = as.numeric(time)
          M     = length(time)
          depth = data[,1]
          data1 = as.matrix(t(data[,2:ncol(data)]))

          if (var == 'Aks'){
             #Calculate MLD based on threshold of Kv (1E-4 m2/s):
             MLD <- numeric(M)
             for (i in 1:M){
                 MLD[i] <- interp(x=depth,y=data[,i+1],Xth=1E-4)
             } 
             #Calculate MLD based on observation profiles:
             MLD_ob <- MLD_obs(Stn_name)
          }
          x     = which(depth>=-250)
          depth = depth[x]

          if (var == 'Aks'){
            data1 = log10(data1[,x])
           title1 = bquote('Log'[10] 
                    ~ italic(' K')[v] ~ ' at ' ~
                    .(Stn_name))
        
          } else{
              data1  <- data1[,x]
              title1 <- var
             if (var=='temp'){
                 title1 <- bquote('Temperature (ºC) at ' ~ .(Stn_name))
             }
          }

          image2D(data1, x=as.numeric(time), y=depth,lwd=2,
                      xlab="Month", ylab="Depth (m)")
          mtext(title1,line = 0.3)
          if (var=='Aks') {
              lines(time,MLD,col='tan',lwd=5)
              points(MLD_ob$DOY/30,MLD_ob$MLD,pch=0,cex=1.5,col='white')
          }
          axis(1, at = seq(1,11,by=2))
       dev.off()
}
