filedir  <- paste('~/working/FlexEFT1D/DRAM_0.2/',
                   Stn,'/',model,sep='')
setwd(filedir)
library(gridExtra)		 
#Plot correlations of the parameters
#Open enspar 
enspar <- read.table('enspar',header=T)
Z      <- enspar[, 3:ncol(enspar)]
nms    <- names(enspar)[3:ncol(enspar)]

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

burnin  = 8E+4
enspar  = enspar[burnin:nrow(enspar),]

 graph <- list()
for (ii in 1:(length(nms)-1)){
   for (jj in (ii+1):length(nms)){

   enspar$x  <- enspar[,nms[ii]]
   enspar$y  <- enspar[,nms[jj]]
   corp      <- cor(enspar$x,enspar$y)
   xx        <- ggplot(enspar, aes(x=x, y=y)) +
          stat_density2d(aes(fill = ..level..), geom="polygon")+
          scale_fill_gradientn(colours=rev(rainbow(100, start = 0.1, end=.8)))+  
          xlab(labs[ii]) +
          ylab(labs[jj]) +
          annotate("text",
           x      = median(enspar$x,na.rm=T), 
           y      = median(enspar$y,na.rm=T), 
           label  = paste('r =',round(corp,2)), 
           family = "serif", size=4)+

          LO_theme 
  # ggExtra::ggMarginal(xx, type = "histogram")
   graph  <- c(graph,  list(xx))
   }
}

ensfile <- paste(Stn,model,'_ens_pars.pdf',sep='')
#ensfile <- paste(Stn,model,'_ens_pars.tiff',sep='')
pdf(ensfile,paper='a4',width=9,height=12)
#tiff(ensfile,width=5,height=6,units='in',res=300)
marrangeGrob(graph,nrow=4,ncol=3)
dev.off()

