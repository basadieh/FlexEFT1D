#Read data files from 1D model output:
setwd('/Users/apple/working/FlexEFT1D/S1/FlexDISC')
library(plot3D)

getData = function(var){
   var  = as.character(var)
   file = paste(var,'.out',sep='')
   data = read.table(file,header=F)
   days = as.numeric(as.character(data[2:nrow(data), 2]))
   depth= as.numeric(data[1,3:ncol(data)])
   data = data[2:nrow(data), 3:ncol(data)]
   return(list(days=days,depth=depth,data=data))
}
# Calculate total Phyto N
Aks = getData('Aks')$data


Par = getData('PAR')$data
NO3 = getData('NO3')
NPHY= 20
PHY   <- array(NA, c(dim(NO3$data),NPHY))
QN    <- PHY
Theta <- PHY
TheHat<- PHY
mu    <- PHY
muHat <- PHY
SI    <- PHY
Lno3  <- PHY

days= NO3$days
depth =NO3$depth
dtdays=300/86400
for (i in 1:NPHY){
  modvar   <- paste('PHY',formatC(i,width=2,flag= ''),sep='')
  PHY[,,i] <- as.matrix(getData(modvar)$data)
  modvar   <- paste('QN',formatC(i,width=3,flag= ''),sep='')
  QN[,,i]  <- as.matrix(getData(modvar)$data)
  modvar   <- paste('The',formatC(i,width=2,flag= ''),sep='')
  Theta[,,i]  <- as.matrix(getData(modvar)$data)

  modvar   <- paste('mu',formatC(i,width=3,flag= ''),sep='')
  mu[,,i]  <- as.matrix(getData(modvar)$data)/dtdays

  modvar   <- paste('THA',formatC(i,width=2,flag= ''),sep='')
  TheHat[,,i]  <- as.matrix(getData(modvar)$data)

  modvar   <- paste('SI',formatC(i,width=3,flag= ''),sep='')
  SI[,,i]  <- as.matrix(getData(modvar)$data)

  modvar   <- paste('LNO',formatC(i,width=2,flag= ''),sep='')
  Lno3[,,i]<- as.matrix(getData(modvar)$data)

}

PMU_min=0.123; PMU_max=10
PMU=seq(PMU_min,PMU_max,length.out=NPHY)
#Plot size distribution
phy1 <- PHY[length(days),length(depth),]
barplot(phy1,names.arg =round(PMU,1),beside=T)

the1<-Theta[length(days),length(depth),]
barplot(the1)

thehat1<-TheHat[length(days),length(depth),]
barplot(thehat1)

mu1 <- mu[length(days),length(depth),]
barplot(mu1,names.arg =round(PMU,1),beside=T)

SI1 <- SI[length(days),length(depth),]
barplot(SI1,names.arg =round(PMU,1),beside=T)

Lno31 <- Lno3[length(days),length(depth),]
barplot(Lno31,names.arg =round(PMU,1),beside=T)

QN1 <- QN[length(days),length(depth),]
barplot(QN1,names.arg =round(PMU,1),beside=T)

Chl  <- PHY/QN*Theta
#Calculate total PHY NITROGEN FOR EACH COMMUNITY
PHYt <- apply(PHY, c(1,2),sum)
Chlt <- apply(Chl, c(1,2),sum)

plot_depth = function(var,Dep){
      if (var %in% c('PMU','VAR')){
      PHY = getData('PHY')$data
      data= getData(var)
      days= data$days
     depth= data$depth
      data= data$data/PHY
   } else{ 
  
      data=getData(var)
      days=data$days
     depth=data$depth
      data=data$data
   }
   Dep = -abs(Dep)
   if (Dep >= max(depth)) {
      x   =data[,ncol(data)]
   }else if(Dep <= min(depth)){
      x   =data[,1]
   }else{
      for (i in 1:length(depth)){
          if ( (Dep >= depth[i]) && (Dep < depth[i+1] )){
             rat=(Dep-depth[i])/(depth[i+1]-depth[i])
             x  =data[,i]*(1-rat)+data[,i+1]*rat
          }  
      }
   }
   plot(days,x, xlab='Days', ylab=var,type='l',xaxt='n')   
   axis(1, at=seq(0,max(days),by=360))
   abline(  v=seq(0,max(days),by=360),lty=3)
}

z=c(-1,-50,-100)
z=-1
for (k in 1:length(z)){

   file=paste('S1_',abs(z[k]),'_m.pdf',sep='')
   pdf(file, width=6,height=9,paper='a4')
   op <- par(font.lab = 1,
              family ="serif",
              mar    = c(4,4,1,3),
              mgp    = c(2,1,0),
              mfrow  = c(4,2))  
   for (i in modvars){
       plot_depth(i,z[k])
   }
   dev.off()
}

plot_1D = function(var){
   if (var %in% c('PMU','VAR')){
      PHY = getData('PHY')$data
      data= getData(var)
      days= data$days
     depth= data$depth
      data= data$data/PHY
   } else if(var == 'PHYt'){
     data= PHYt
   } else if(var == 'Chlt'){
     data= Chlt
   } else{ 
      data=getData(var)
      days=data$days
     depth=data$depth
      data=data$data
   }
    image2D(as.matrix(data), x=days, y=depth, 
           xlab='Days',ylab='Depth (m)',main=var,xaxt='n')
   axis(1, at=seq(0,max(days),by=360))
   abline(  v=seq(0,max(days),by=360),lty=3)
}
plot_1D('NO3')
plot_1D('PHYt')
plot_1D('Chlt')

modvars  = c( 'Temp ','PAR  ','Aks  ','NO3  ','PHY  ','ZOO  ','DET  ','PMU  ',
      'VAR  ','CHL  ','Theta','QN   ','muNet','Graz ','Z2N  ',
      'D2N  ','dmudl','d2mu ','d2gdl',
      'D_NO3','D_PHY','D_ZOO','D_DET','D_PMU','D_VAR','D_CHL') 

modvars = gsub(" ", "", modvars, fixed = TRUE)

pdf('S1.pdf', width=6,height=9,paper='a4')
op <- par(font.lab = 1,
           family ="serif",
           mar    = c(4,4,1,3),
           mgp    = c(2,1,0),
           mfrow  = c(4,2))  
for (i in modvars){
    plot_1D(i)
}
dev.off()





NO3 = NO3$data
ZOO = getData('ZOO')$data
DET = getData('DET')$data

Hz <- function(hmax=250, thetaS=2, nlev=40){
  Z_w=numeric(nlev+1)
  Z_r=numeric(nlev)
  Hz =numeric(nlev)
  Z_w[1] = -hmax
  for( i   in 1:nlev){
      sc_r    = ((i-nlev) - 0.5)/nlev
     C_sig    = sinh(thetaS*sc_r)/sinh(thetaS)      # -1 < C_sig < 0
     Z_r[i]   = C_sig*hmax
      sc_r    = (i-nlev)/nlev
     C_sig    = sinh(thetaS*sc_r)/sinh(thetaS)      # -1 < C_sig < 0
     Z_w[i+1] = C_sig*hmax
      Hz[i]   = Z_w[i+1] - Z_w[i]
  }
  return(Hz)
}
Hz <- Hz()


#Calculate total N for each date:
TNO3 = numeric(length(days)) 
TZOO = numeric(length(days)) 
TDET = numeric(length(days)) 
TPHY = numeric(length(days)) 
for (i in 1:length(days)){

     
    TNO3[i] <- sum(NO3[i,]*Hz)   #Unit: mmol/m2
    TDET[i] <- sum(DET[i,]*Hz)   #Unit: mmol/m2
    TZOO[i] <- sum(ZOO[i,]*Hz)   #Unit: mmol/m2
    TPHY[i] <- sum(PHYt[i,]*Hz)   #Unit: mmol/m2

}

#Calculate NO3 budget:

#1) Calculate phyto uptake = PP(iPHY,iNO3):
PHY = getData('PHY')$data
mu  = getData('muNet')$data[,N]

#Calculate surface:
N   = dim(PHY)[2]
N2P = PHY[,N]*mu

#2) Calculate zooplankton excretion:
Z2N = getData('Z2N')$data[,N]
 
#3) Calculate detritus regeneration:
D2N  = getData('D2N')$data[,N]

#4) Calculate NO3 diffusion:
D_NO3=getData('D_NO3')$data[,N]


plot(days,N2P,type = 'l',xaxt='n')
axis(1, at=seq(0,max(days),by=360))
abline(  v=seq(0,max(days),by=360),lty=3)

points(days,Z2N,type='l',col=2)
points(days,D2N,type='l',col='green')
points(days,D_NO3,type='l',col='blue')

dN   <- (Z2N + D2N + D_NO3) - N2P
NO3  <- getData('NO3')$data[N]
NO3  <- as.vector(NO3)
dN2  <- NO3[2:nrow(NO3),] - NO3[1:(nrow(NO3)-1),]
dN2  <- c(0,dN2) 
plot(days,dN,type='l')


plot(dN,dN2,xlim=c(-0.005,0.005),ylim=c(-0.005,0.005))

