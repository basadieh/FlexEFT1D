#Plot LogLikehood of two models at two stations in one plot:
setwd('~/Working/FlexEFT1D/')  
#pdffile <- paste('LogLike.pdf',sep='')
#pdf(pdffile, width=5,height=5,paper='a4')
tiff('~/Working/FlexEFT1D/LogLike.tiff',width=5,height=5,units='in',res=300)
op <- par(font.lab = 1,
           family ="serif",
           mar    = c(2,4,2,.3),
           mgp    = c(2,1,0),
           oma    = c(2.3,0,0,0),
           mfrow  = c(4,4))  
ii=0
for (stn in c('S1','K2')){
  for (model in c('Geidersimple','EFTsimple')){ 
      burnin  = 1E+4
      NDTYPE  = 3
      filedir = paste('~/working/FlexEFT1D/DRAM_0.2/',stn,'/',model,sep='')
      setwd(filedir)
      enssig  <- read.table('enssig',header=T)
      enssig  <- enssig[burnin:nrow(enssig),]
      runno   <- enssig[,1]
      loglike <- enssig[,2]
      sigma   <- enssig[,(2+1):(2+NDTYPE)]
      ssqe    <- enssig[,(ncol(enssig)-NDTYPE+1):ncol(enssig)]
      
      plot(runno/1E+4,loglike,
            xlab='',
            ylab='Log-likelihood',
            cex =0.1,pch=16)

      ii=ii+1
      mtext(letters[ii],adj=0,cex=.8)
      model_ = rep(c(rep('NPZDChl',4), rep('NPZDFlex',4)),2)
      if (ii%%4 == 1) mtext(paste(stn,model_[ii]),adj=.6,cex=.6)

      varnames <- c( 'TIN','Chl','NPP')
      for (i in 1:ncol(ssqe)){
          plot(runno/1E+4,ssqe[,i],
               xlab='',
               ylab=paste('SSqE_',varnames[i],sep=''),
               cex =0.1,pch=16)
          ii=ii+1
          mtext(letters[ii],adj=0,cex=1)
      }
  } 
}
xlab  = expression(paste('Run number ('*10^4*')',sep=''))
mtext(xlab,side=1,outer=TRUE,cex=.8)
dev.off()

