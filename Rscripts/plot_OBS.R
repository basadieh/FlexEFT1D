#Plot Chl
setwd('~/Working/FlexEFT1D')
library(mgcv)
library(plot3D)
source('Rscripts/surface_Chl.R')
plot_OBS <- function(var,Stn,method='linear'){
   setwd(paste('~/Working/FlexEFT1D/',Stn,sep=''))  
   dat   = paste(Stn,'_',var,'.dat',sep='')
   if (file.exists(dat)){
     dat = read.table(dat,header=T)
   } else{
     stop(paste(dat,'does not exist!'))
   }
   dat$y =dat[,3]
   dat$y[dat$y <= 0] = .001
   range =quantile(dat$y,probs=c(0.05,0.95))
   range =round(range,1)
   #Construct new matrix
   newD  =seq(250,1,-1)
   newDOY=seq(1,360,1)
  
   if (method == 'GAM'){
      newMat=expand.grid(DOY=newDOY,Depth=newD)
      cff   =gam(log(y) ~ te(DOY,Depth), data=dat, gamma=1.4)
      cff   =predict.gam(cff,newdata=newMat)
      cff   =exp(cff)
      w     =which(newMat$Depth > max(dat$Depth))
      if(var %in% c('Chl','NPP')){
         cff[w]=.001
      }
      z       <- data.frame(matrix(cff, nc=length(newDOY), nr=length(newD),byrow=T))
      pdfname <- paste(Stn,'_',var,'_gam.pdf',sep='')
     }else if(method == 'linear'){
      #Remove only surface data:
      zz =tapply(dat$Depth,dat$DOY,length)
      zz1=tapply(dat$Depth,dat$DOY,max)
 
      #Determine unique profiles:
      DOY_uni = unique(dat$DOY)
      DOY_uni = sort(DOY_uni)
      DOY_uni = DOY_uni[zz>3 & zz1>100]
      #Construct an intermediate matrix to store the data in newD
      dmy     = data.frame(matrix(NA,nr=length(newD),nc=length(DOY_uni)))

      for (i in 1:length(DOY_uni)){
          DOY_ = DOY_uni[i]
          pf   = dat[dat$DOY == DOY_,]
          pf   = pf[order(pf$Depth),]

          #For each profile, linearly interpolate to the targeted depth in newD
          Dep_min = min(pf$Depth)
          Dep_max = max(pf$Depth)

          #Unique profile depths:
          pf_deps = unique(pf$Depth)

          if (Dep_max > 20){
             for (k in 1:length(newD)){
                 #Targeted depth
                 Dep_t = newD[k]
                 if (Dep_t >= Dep_max) {
                    w        = which(pf$Depth == max(pf$Depth,na.rm=T))
                    dmy[k,i] = mean(pf[w,ncol(pf)],na.rm=T)
                 } else if(Dep_t <= Dep_min){
                    w        = which(pf$Depth == min(pf$Depth,na.rm=T))
                    dmy[k,i] = mean(pf[w,ncol(pf)],na.rm=T)
                 } else{
                    for (j in 1:(length(pf_deps)-1)){
                        pf_dep1 = pf_deps[j  ]
                        pf_dep2 = pf_deps[j+1]
                        if (pf_dep1 < Dep_t && pf_dep2 >= Dep_t){
                            #Find the corresponding y values in pf:
                            y1 = pf[pf$Depth == pf_dep1,ncol(pf)]
                            y1 = mean(y1,na.rm=T)
                            y2 = pf[pf$Depth == pf_dep2,ncol(pf)]
                            y2 = mean(y2,na.rm=T)
                            ra = (Dep_t-pf_dep1)/(pf_dep2-pf_dep1)
                            dmy[k,i] = y1*(1-ra)+y2*ra
                            break
                        }

                    }
                 }

             }
          }
      }

      #For each depth,   linearly interpolate dmy to the targeted DOY   in newDOY
       z  <- data.frame(matrix(NA, nc=length(newDOY), nr=length(newD),byrow=T))

         w_min <- which.min(DOY_uni)
         w_max <- which.max(DOY_uni)
       DOY_min <- min(DOY_uni)
       DOY_max <- min(c(max(DOY_uni),360))

       # The data corresponding to minimal DOY:
       zmin <- dmy[,w_min]
       zmin <- matrix(zmin,nc=1)
       # The data corresponding to maximal DOY:
       zmax <- dmy[,w_max]
       zmax <- matrix(zmax,nc=1)

             
       #Using loess to interpolate and smooth the data:
       for (k in 1:length(newD)){
           #Obs data for the targeted depth:
           y     <- dmy[k,]
           y     <- data.frame(DOY=DOY_uni,dat=as.numeric(y))
           fit   <- loess(dat ~ DOY, y)
           newy  <- predict(fit, data.frame(DOY = newDOY))
           z[k,] <- as.numeric(newy)
       }

       #Interpolate data at both bounds
       for (j in 1:ncol(z)){
           #Targeted DOY:
           DOY_t   <- newDOY[j]
           if (DOY_t <= DOY_min){
              ra     <- (DOY_t-DOY_max+360)/(DOY_min-DOY_max+360)
              z[,j]  <- zmin*ra + zmax*(1-ra)
           }else if(DOY_t >= DOY_max){
              ra     <- (DOY_t-DOY_max)/(DOY_min-DOY_max+360)
              z[,j]  <- zmin*ra + zmax*(1-ra)
           }else{
              # DO NOTHING
              # for (ii in 1:(length(DOY_uni)-1)){
              #     if (DOY_t >= DOY_uni[ii] && DOY_t < DOY_uni[ii+1]){
              #        ra     <- (DOY_t-DOY_uni[ii])/(DOY_uni[ii+1]-DOY_uni[ii])
              #        z[,j]  <- dmy[,ii]*(1-ra)+dmy[,ii+1]*ra
              #        break
              #     }
              # }
           } 
       }
       


      pdfname <- paste(Stn,'_',var,'_linear.pdf',sep='')
   }  # of if method == 'linear'
   z             <- t(z)
   z[z>range[2]] <- range[2]
   z[z<range[1]] <- range[1]

   pdf(file=pdfname,height=5,width=5,paper='a4')
       par(font.lab  = 1,
           family    = "serif",
           mgp       = c(2.2,1,0))

       image2D(z,x=newDOY,y=newD,zlim=range,lwd=2,
               xlab='Date of the year',ylab='Depth (m)',xaxt='n')
       points(dat$DOY,dat$Depth,pch=16,cex=.5)
       axis(1, at=seq(0,360,by=60))
       title(main=paste('Observed',var,'at',Stn))
   dev.off()
   setwd('~/Working/FlexEFT1D')
}


plot_OBS('TIN','S1')
plot_OBS('Chl','S1')
plot_OBS('NPP','S1')
plot_OBS('TIN','K2')
plot_OBS('Chl','K2')
plot_OBS('NPP','K2')
