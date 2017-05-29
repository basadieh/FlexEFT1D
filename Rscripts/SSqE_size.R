#Plot SSqE of two models at two stations in one plot:
setwd('~/Working/FlexEFT1D/')  
pdffile <- paste('ssqe_size.pdf',sep='')
pdf(pdffile, width=5,height=5,paper='a4')
op <- par(font.lab = 1,
           family ="serif",
           mar    = c(2,2,2,.3),
           mgp    = c(2,1,0),
           oma    = c(0,2,0,0),
           mfrow  = c(2,2))  
ii=0
for (stn in c('S1','K2')){
  for (model in c('EFTdiscrete','EFTcont')){ 
   ii=ii+1
   burnin  = 8E+3
   NDTYPE  = 3
   Ymin    = 9
   Ymax    = 35
   filedir = paste('~/working/FlexEFT1D/DRAM_0.2/',stn,'/',model,sep='')
   setwd(filedir)
   enssig  <- read.table('enssig',header=T)
   enssig  <- enssig[burnin:nrow(enssig),]
   ssqe    <- enssig[,(ncol(enssig)-NDTYPE+1):ncol(enssig)]
   
   model_  <- rep(c('NPZDChl', 'NPZDFlex'),2)
   boxplot(ssqe,names=c('TIN','Chl','NPP'),
           ylim    = c(Ymin, Ymax),
           cex.lab = 1.2, cex.axis = 1.2,
           ylab    = '',
#           ylab   = 'Squared errors of transformed variables',
           outline = F,col = "lightgray")  
   mtext(letters[ii],adj=0,cex=1)
   mtext(paste(stn,model_[ii]),adj=.6,cex=1)
  } 
}
ylab   = 'Squared errors of transformed variables'
mtext(ylab,side=2,outer=TRUE)
dev.off()

