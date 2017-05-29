dat <- read.csv2('S1_obs.csv')
dat$Lon   = as.numeric(as.character(dat$Longitude))
dat$Lat   = as.numeric(as.character(dat$Latitude))
dat$Date  = as.Date(dat$Date,'%m/%d/%Y')
dat$DOY   = as.numeric(strftime(dat$Date, format = "%j"))    #Date of the year
dat$Depth = as.numeric(as.character(dat$Depth)) 

#Clean up nitrate data:
dat$Nitrate=as.numeric(as.character(dat$Nitrate))
dat$Nitrate[dat$Nitrate<0]=NA
dat$Nitrate1=as.numeric(as.character(dat$Nitrate1))
dat$Nitrate1[dat$Nitrate1<0]=NA

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


dat$NO3   = sapply(1:nrow(dat),function(i)mean(c(dat$Nitrate[i], dat$Nitrate1[i]),na.rm=T)) 
dat$NO2   = sapply(1:nrow(dat),function(i)mean(c(dat$Nitrite[i], dat$Nitrite1[i]),na.rm=T)) 
dat$Nh4   = sapply(1:nrow(dat),function(i)mean(c(dat$NH4[i], dat$NH41[i]),na.rm=T)) 

#Add up
dat$TIN   = sapply(1:nrow(dat),function(i)sum(c(dat$NO3[i], dat$NO2[i], dat$Nh4[i]),na.rm=T)) 

#Clean up CHL data
dat$CHL=as.numeric(as.character(dat$CHLWEL))
dat$CHL[dat$CHL<0]=NA
dat$CHL1=as.numeric(as.character(dat$CHLWEL1))
dat$CHL1[dat$CHL1<0]=NA
dat$Chl   = sapply(1:nrow(dat),function(i)mean(c(dat$CHL[i], dat$CHL1[i]),na.rm=T)) 
dat$Chl[(dat$Chl <= 0)& (!is.na(dat$Chl))] = 0.01 

#Clean up TIN data

dat$TIN[dat$TIN<0]=NA
dat$TIN[(dat$TIN <= 0)& (!is.na(dat$TIN))] = 0.01 

#Only need data above 250 m
dat = dat[dat$Depth <= 250,]

#Plot Chl
library(mgcv)
library(plot3D)

plot_OBS <- function(var,Stn){
   dat$y = dat[,var]
   gam2D = gam(log(y) ~ te(Depth, DOY),data=dat,gamma=1.4)
   
   #Construct new matrix
   newD  =seq(250,1,-1)
   newDOY=seq(1,360,1)
   newMat=expand.grid(Depth=newD,DOY=newDOY)
   newMat$y=exp(predict.gam(gam2D,newdata=newMat))
   z  <- matrix(newMat$y, nc=length(newDOY), nr=length(newD))
   z  <- t(z)
   pdfname <- paste(Stn,'_',var,'.pdf',sep='')
   pdf(file=pdfname,height=5,width=5,paper='a4')
       par(font.lab  = 1,
           family    = "serif",
           mgp       = c(2.2,1,0))

       image2D(z,x=newDOY,y=newD,xlab='Date of the year',ylab='Depth (m)')
       points(dat$DOY,dat$Depth,pch=16,cex=.8)
       title(main=paste('Observed',var,'at',Stn))
   dev.off()
}


plot_OBS('TIN','S1')
plot_OBS('Chl','S1')
