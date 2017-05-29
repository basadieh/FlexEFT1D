#Plot the model output of surface Chl:
source('~/Working/FlexEFT1D/Rscripts/getData.R')

pdf('Surface_mod_obs.pdf', width=5,height=6,paper='a4')
   op <- par(font.lab = 1,
                family ="serif",
                mar    = c(4,4,1,3),
                mgp    = c(2,1,0),
                mfrow  = c(1,1)) 

for (Stn in c('S1','K2')){
  for (model in c('EFTsimple','Geidersimple')){
     filedir  <- paste('~/working/FlexEFT1D/DRAM_0.2/',
                        Stn,'/',model,sep='')
     setwd(filedir)

     #Read observational data
     dat     = paste('~/Working/FlexEFT1D/',Stn,'/',Stn,'_CHL.dat',sep='')
     dat     = read.table(dat,header=T)
     dat     = dat[dat$Depth <= 10,]
     ymax    = quantile(dat$Chl,probs=0.99)
     #Plot data points:
     plot(dat$DOY,dat$Chl,
       xlab='Date of the year',
       ylab=expression(paste("Chl "*italic(a)*' (µg '*L^-1*')')),
       xlim=c(0,360),ylim=c(0,ymax),cex=.3,pch=16,xaxt='n')

     #Get modeled Chl data
        Chl  = getData('CHL_T')
        days = Chl$days
       depth = Chl$depth
        Chl  = Chl$data
        d_per_y = 360
        w       = (length(days)-d_per_y+1):length(days)
        Chl     = Chl[w,]
        w       = depth >= -10
        s_Chl   = Chl[,w]
        s_Chl   = apply(s_Chl,1,mean)
    
        lines(1:length(s_Chl),s_Chl, lwd=1.5)


  }
}
