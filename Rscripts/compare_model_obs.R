
setwd('~/Working/FlexEFT1D/')  
ii=0
pdffile <- paste('Model_obs_compare.pdf',sep='')
pdf(pdffile, width=5,height=6,paper='a4')
op <- par(font.lab = 1,
           family ="serif",
           mar    = c(3,3,1,.3),
           mgp    = c(2,1,0),
           oma    = c(0,0,0,0),
           mfrow  = c(4,3))  

for (stn in c('S1','K2')){
  for (model in c('Geidersimple','EFTsimple')){ 

   filedir  <- paste('~/working/FlexEFT1D/DRAM_0.2/',
                       stn,'/',model,sep='')
   setwd(filedir)
   # Compare bestoutput with data:
   bestout <- read.table('bestout',header=T)
   for (var in c('TIN','CHL','NPP')){
       ii=ii+1
       obsfile <- paste(stn,'_',var,'.dat',sep='')
       obs     <- read.table(obsfile,header=T)
       mod     <- bestout[bestout$Data_type == var,]
       x       <- obs[,ncol(obs)]
       y       <- mod[,ncol(mod)]
       Xmax    <- as.numeric(quantile(x,probs=0.95))
       model_  <- rep(c(rep('NPZDChl',3), rep('NPZDFlex',3)),2)
       plot(x,y,
            xlim=c(0,Xmax),
            ylim=c(0,Xmax),
            xlab=paste('Observed',var),
            ylab=paste('Modeled',var),
            cex =.2, pch=1)
       abline(0,1,lty=2,lwd=.5)
       text(0,Xmax*.95,cex=1,pos=4,paste('r =',round(cor(x,y),2)))
       mtext(letters[ii],adj=0,cex=1)
       if (ii%%3 == 1) mtext(paste(stn,model_[ii]),adj=.6,cex=.8)
   }
  }
}
dev.off()

