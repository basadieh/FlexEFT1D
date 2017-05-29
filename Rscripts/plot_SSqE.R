
source('~/Working/FlexEFT1D/Rscripts/param_lab.R')
plot_SSqE <- function(Stn,model_ID,burnin=100){
     filedir  <- paste('~/working/FlexEFT1D/DRAM_0.2/',
                        Stn,'/',model_ID,sep='')
     setwd(filedir)

     if (model_ID == 'Geidersimple') {
      labs  <- c(mu0,alphaI,Q0N,KN,Gm,Kp,Wdet,Rdn,mz)
      NDTYPE<- 3
     } else if(model_ID == 'EFTsimple'){
      labs  <- c(alphaI,A0N,Q0N,Gm,Kp,Wdet,Rdn,mz)
      NDTYPE<- 3
     } else if(model_ID == 'EFTsimple1'){
      labs  <- c(alphaI,KN,Q0N,Gm,Kp,Wdet,Rdn,mz)
      NDTYPE<- 3
     }

     enspar  <- read.table('enspar',header=T)
     enspar  <- enspar[burnin:nrow(enspar),]
     Z       <- enspar[,3:ncol(enspar)]
     nms     <- names(enspar)[3:ncol(enspar)]
     #open enssigma:
     enssig  <- read.table('enssig',header=T)
     enssig  <- enssig[burnin:nrow(enssig),]
     runno   <- enssig[,1]
     loglike <- enssig[,2]
     sigma   <- enssig[,(2+1):(2+NDTYPE)]
     ssqe    <- enssig[,(ncol(enssig)-NDTYPE+1):ncol(enssig)]

     #Check if convergency has been reached:
     pdffile <- paste(Stn,model_ID,'_1D_param.pdf',sep='')
     pdf(pdffile, width=8,height=8,paper='a4')
     op <- par(font.lab = 1,
                family ="serif",
                mar    = c(4,4,1,3),
                mgp    = c(2,1,0),
                mfrow  = c(2,2)) 
     plot(runno,loglike,
              xlab='Run number',
              ylab='Log-likelihood',
              cex =0.4,pch=16)

     varnames <- c( 'TIN','Chl','NPP')
     for (i in 1:ncol(ssqe)){
         plot(runno,ssqe[,i],
              xlab='Run number',
              ylab=paste('SSqE_',varnames[i],sep=''),
              cex =0.4,pch=16)
     }
     for (i in 1:ncol(sigma)){
         plot(runno,sigma[,i],
              xlab='Run number',
              ylab=paste('Sigma_',varnames[i],sep=''),
              cex =0.4,pch=16)
     }

     for (i in 1:length(nms)){
         plot(enspar$Run,Z[,i],
              xlab='Run number',
              ylab=labs[i],cex=0.4,pch=16) 
     }
     dev.off()

     #write average (and SD) parameter values into a file 
     mat <- data.frame(matrix(NA, nrow = 2, ncol = length(nms)))
     names(mat) <- nms
     mat[1,]    <- colMeans(Z)
     mat[2,]    <- sapply(1:ncol(mat), function(i)sd(Z[,i]))
     write.table(round(mat,3), file=paste(Stn,model_ID,'_meanSDpar',sep=''),
                append=F,row.names=F,col.names=T)  
     #Plot SSqEs:
     pdffile <- paste(Stn,model_ID,'_ssqe.pdf',sep='')
     pdf(pdffile, width=5,height=5,paper='a4')
     op <- par(font.lab = 1,
                family ="serif",
                mar    = c(4,4,1,3),
                mgp    = c(2,1,0),
                mfrow  = c(1,1))  

     boxplot(ssqe,names=c('TIN','Chl','NPP'),
             cex.lab = 1.2, cex.axis = 1.2,
             ylab    = 'Squared errors of transformed variables',
             outline = F,col = "lightgray")  
     dev.off()
}

