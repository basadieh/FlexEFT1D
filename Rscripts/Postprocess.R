#Read data files from 1D model output:
library(plot3D)
source('~/Working/FlexEFT1D/Rscripts/LO_theme.R')
source('~/Working/FlexEFT1D/Rscripts/param_lab.R')
#Compare SSqE among model and station
source('~/Working/FlexEFT1D/Rscripts/SSqE.R')
#Plot changes of LogLikelihood and SSqE
source('~/Working/FlexEFT1D/Rscripts/loglike.R')
source('~/Working/FlexEFT1D/Rscripts/loglike_size.R')  #Plot Log-Likelihoods of size model
source('~/Working/FlexEFT1D/Rscripts/compare_model_obs.R')
source('~/Working/FlexEFT1D/Rscripts/compare_param.R')

#++++++++++++++++++++++Dynamics within the MLD+++++++++++++++++++++++
#plot out data of CHL and NO3 and model fits in the mixed layer
source('~/Working/FlexEFT1D/Rscripts/surface_OBS_model.R')
source('~/Working/FlexEFT1D/Rscripts/MLD_OBS_model_size.R')   #For the size models

#Plot out data of PHY, ZOO and Detritus in the mixed layer
source('~/Working/FlexEFT1D/Rscripts/MLD_PHYZOO.R')

#plot C:Chl ratio and growth/grazing rates within the mixed layer
source('~/Working/FlexEFT1D/Rscripts/CChl_ML.R')

#plot responses of growth rates and C:Chl ratio to light and nutrient for two models
source('~/Working/FlexEFT1D/Rscripts/mu_theta_response.R')

#For two size models:
#plot mean size and variance at the two stations
source('~/Working/FlexEFT1D/Rscripts/MLD_size_var.R')

#plot phytoplankton size histograms at two stations:
source('~/Working/FlexEFT1D/Rscripts/Phy_hist.R')
#+++++++++++++++++++++End of MLD scripts+++++++++++++++++++++++++++++

#++++++++++++++++++++Vertical distributions++++++++++++++++++++++++++
#Plot vertical distributions of Chl, NO3, and NPP
#Separation into 4 seasons are needed to reduce the complexity
source('~/Working/FlexEFT1D/Rscripts/vertical_OBS_model.R')


#+++++++++++++++++++++End of Vertical distributions++++++++++++++++++


dtdays  <- 300/86400


Hz <- Hz()

Stns     <- c('S1','K2')
modelIDs <- c('EFTsimple','Geidersimple')

for (Stn in c('S1','K2')){
  for (model in c('EFTsimple','Geidersimple')){
    source('~/Working/FlexEFT1D/Rscripts/param_distr.R')
  }
}

source('~/Working/FlexEFT1D/Rscripts/plot_SSqE.R')
plot_SSqE('S1','Geidersimple')
     

plot_1D('CHL_T')
NPHY     <- 1
for (i in 1:NPHY){
     phyvar   <- paste('PHY',formatC(1,width=2,flag= ''),sep='')
     Thetavar <- paste('The',formatC(1,width=2,flag= ''),sep='')
     muvar    <- paste('mu', formatC(1,width=3,flag= ''),sep='')
     QNvar    <- paste('QN', formatC(1,width=3,flag= ''),sep='')

     plot_1D(phyvar)
     plot_1D(muvar)
     plot_1D(QNvar)
     plot_1D(Thetavar)

}


var     <- 'TIN'
obsfile <- paste(Stn,'_',var,'.dat',sep='')
obs     <- read.table(obsfile,header=T)
mod     <- bestout[bestout$Data_type == var,]
x       <- obs[,ncol(obs)]
y       <- mod[,ncol(mod)]

obsmin <- min(x)
obsmax <- max(x)
x      <- transf(x,obsmin,obsmax)
y      <- transf(y,obsmin,obsmax)
hist(x)


# Calculate total Phyto N
Aks = getData('Aks')$data


Par = getData('PAR')$data
PHY   <- array(NA, c(dim(NO3$data),NPHY))
QN    <- PHY
Theta <- PHY
TheHat<- PHY
mu    <- PHY
muHat <- PHY
SI    <- PHY
Lno3  <- PHY

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

plot_1D('PHYt')
plot_1D('Chlt')

modvars  = c( 'Temp ','PAR  ','Aks  ','NO3  ','PHY  ','ZOO  ','DET  ','PMU  ',
      'VAR  ','CHL  ','Theta','QN   ','muNet','Graz ','Z2N  ',
      'D2N  ','dmudl','d2mu ','d2gdl',
      'D_NO3','D_PHY','D_ZOO','D_DET','D_PMU','D_VAR','D_CHL') 

modvars = gsub(" ", "", modvars, fixed = TRUE)

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

