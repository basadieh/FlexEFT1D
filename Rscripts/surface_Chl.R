library(plot3D)
library(mgcv)
#plot seasonal changes of surface Chl

plot_surChl = function(Stn){
   setwd(paste('~/Working/FlexEFT1D/',Stn,sep=''))  
   dat   = paste(Stn,'_CHL.dat',sep='')
   if (file.exists(dat)){
     dat = read.table(dat,header=T)
   } else{
     stop(paste(dat,'does not exist!'))
   }
   dat=dat[dat$Depth <= 10,]
   z  =gam(log(Chl) ~ s(DOY,bs='cc'), data=dat, gamma=1.4)
   newx=data.frame(DOY=1:360)
   newy=exp(as.numeric(predict(z,newx)))
   
   ymax=quantile(dat$Chl,probs=0.95)
   pdffile <- paste('Surface_Chl_',Stn,'.pdf',sep='')
   pdf(pdffile,width=3.5,height=4,paper='a4')
   par(font.lab  = 1,
       family    = "serif",
       mgp       = c(2.2,1,0),
       mfrow     = c(1,1))

   plot(dat$DOY,dat$Chl,
        xlab='Date of the year',
        ylab=expression(paste("Chl "*italic(a)*' (µg/L)')),
        xlim=c(0,360),ylim=c(0,ymax),cex=.3,pch=16,xaxt='n')
   lines(newx$DOY,newy,lwd=2)
   axis(1, at=seq(0,360,by=60))
   dev.off()
   setwd('~/Working/FlexEFT1D')
}
