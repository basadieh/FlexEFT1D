plot_1D <- function(var,title,ZLIM=NULL, finalyr = F){
 Theta     <- paste('The', formatC(1,width=2,flag= ''),sep='')
 file      <- paste(var,'.out',sep='')

 if (!file.exists(file)) {
   if (var == Theta){
      #Calculate Chl:C ratio from Geidersimple model
      PHY  <- getData('PHY 1')
      days <- PHY$days
     depth <- PHY$depth
      PHY  <- PHY$data
     #Get Chl data
      Chl  <- getData('CHL_T')$data
     #Get Q0N (N:C ratio) from bestpar
      bestpar <- read.table('bestpar',skip=5)
      kk      <- which(bestpar[,1] == 'Q0N')
      Q0N     <- bestpar[kk,4]
     #Calculate Theta:
      data <- Chl/(PHY/Q0N)  #Unit: gChl/(molC)
   } else if(var == 'LNO 1'){
     #Calculate N limitation index from Geidersimple model 
     #Get NO3 data
      NO3     <- getData('NO3')
       days   <- NO3$days
     depth    <- NO3$depth
      NO3     <- NO3$data

     #Get KN from bestpar
      bestpar <- read.table('bestpar',skip=5)
      kk      <- which(bestpar[,1] == 'KN')
      KN      <- bestpar[kk,4]
      data    <- 1-NO3/(NO3 + KN) 

   } else if(var == 'SI  1'){
     #Calculate light limitation from Geidersimple model
     #Get aI0 from bestpar
      bestpar <- read.table('bestpar',skip=5)
      kk      <- which(bestpar[,1] == 'aI0')
      aI0     <- bestpar[kk,4]
     #Get PAR:
      PAR     <- getData('PAR')
      days    <- PAR$days
     depth    <- PAR$depth
      PAR     <- PAR$data
     #calculate theta
      PHY  <- getData('PHY 1')
      PHY  <- PHY$data
     #Get Chl data
      Chl  <- getData('CHL_T')$data
     #Get Q0N (N:C ratio) from bestpar
      kk   <- which(bestpar[,1] == 'Q0N')
      Q0N  <- bestpar[kk,4]
     #Calculate Theta:
      theta <- Chl/(PHY/Q0N)  #Unit: gChl/(molC)
     #Calculate rmax_T
      #Get temperature
      temp <- getData('Temp')$data
      #Get mu0 from bestpar
      kk   <- which(bestpar[,1] == 'mu0hat')
      mu0  <- bestpar[kk,4]
      TEMPBOL <- function(tC,Ep = 0.4,kb= 8.62E-5,Tr=15){
        return( exp(-(Ep/kb)*(1/(273.15 + tC)-1/(273.15 + Tr))))
      }
    rmax_T <- mu0* TEMPBOL(temp)
      data <- 1 - exp(-aI0 *PAR * theta/rmax_T)
   }
 } else{
      data  <- getData(var)
      days  <- data$days
     depth  <- data$depth
      data  <- data$data
      if (var == 'LNO 1') data <- data*2 
 }
 
     if (finalyr) {
        #Take the final year   
      fday <- days[length(days)]  #Final day
      cff  <- which( (days > fday-360) & (days <= fday))
      days <- days[cff]%%360
      days[days==0] <- 360
      data <- data[cff,]
    }
    if (!is.null(ZLIM)){
       data[data < ZLIM[1]] <- ZLIM[1]
       data[data > ZLIM[2]] <- ZLIM[2]
       
    }
    image2D(as.matrix(data), x=days, y=-depth,zlim=ZLIM,
           xlab='Days',ylab='Depth (m)',main=title,xaxt='n')
    #Plot nutricline
    if(var == 'NO3') {
      NO3cln <- numeric(nrow(data))
      for (i in 1:nrow(data)){
          x  <- data[i,]
          z  <- which.min(x>1)
      NO3cln[i] <- depth[z]
      }
      lines(days,-NO3cln,lwd=2,col='green')
    }

   if (!finalyr) {
      axis(1, at=seq(0,max(days),by=360))
      abline(  v=seq(0,max(days),by=360),lty=3)
   }else {
      axis(1, at=seq(15,360,by=90), labels=c('Jan','Apr','Jul','Oct'))
   }
}


