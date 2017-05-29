#Compare parameters of the same model at the two stations:
setwd('~/Working/FlexEFT1D/')  
source('Rscripts/param_lab.R')
for (model in c('Geidersimple','EFTsimple')){ 
  if (model == 'Geidersimple') {
   labs  <- c(mu0,alphaI,Q0N,KN,Gm,Kp,Wdet,Rdn,mz)
   NDTYPE<- 3
  } else if(model == 'EFTsimple'){
   labs  <- c(alphaI,A0N,Q0N,Gm,Kp,Wdet,Rdn,mz)
   NDTYPE<- 3
  } else if(model == 'EFTsimple1'){
   labs  <- c(alphaI,KN,Q0N,Gm,Kp,Wdet,Rdn,mz)
   NDTYPE<- 3
  }
 NPar   = length(labs)
 params = array(NA, dim=c(NPar,2,2))

 pdffile <- paste('~/Working/FlexEFT1D/Params_',model,'.pdf',sep='')
 pdf(pdffile, width=5,height=6,paper='a4')
 op <- par(font.lab = 1,
            family ="serif",
            mar    = c(3,4,1,.3),
            mgp    = c(2.4,1,0),
            oma    = c(0,0,0,0),
            mfrow  = c(3,3))  
 
  ii=0
  for (stn in c('S1','K2')){
   ii=ii+1

   filedir  <- paste('~/working/FlexEFT1D/DRAM_0.2/',
                       stn,'/',model,sep='')
   setwd(filedir)
   enspar  = read.table('enspar',header=T)
   burnin  = 8E+4
   enspar  = enspar[burnin:nrow(enspar),]
   Z       = enspar[, 3:ncol(enspar)]
   nms     = names(enspar)[3:ncol(enspar)]
   stopifnot(length(nms)==NPar)
   for (i in 1:length(nms)){
       params[i,1,ii] = mean(Z[,i])
       params[i,2,ii] =   sd(Z[,i])
   }
  }
  ii=0
  for (i in 1:NPar){
      ii=ii+1
      Y1=params[i,1,]-2*params[i,2,]
      Y2=params[i,1,]+2*params[i,2,]
      barcenters=barplot(params[i,1,],
                 names.arg=c('S1','K2'),
                 ylim     =c(0,max(Y2)),
                 ylab     =labs[i],col='white',
                 las=1,beside = true,border = "black", axes = TRUE)
      segments(x0=barcenters,y0=Y1,x1=barcenters,y1=Y2)
      mtext(letters[ii],adj=0,cex=.8)
      box()
  }
  dev.off()
}


