#Compare C:Chl ratios and growth/grazing rates between the two models within ML
setwd('~/Working/FlexEFT1D/')  
source('Rscripts/get_MLD.R')
source('~/Working/FlexEFT1D/Rscripts/interp.R')
source('~/Working/FlexEFT1D/Rscripts/getData.R')
source('~/Working/FlexEFT1D/Rscripts/Interpolate_WOA.R')
#Estimate C:Chl from remote sensing:
Chl_C = readnc('Chl_C')  #A global Chl:C data on the surface
source('~/Working/FlexEFT1D/Rscripts/getdata_station.R')

get_theta = function(Stn_name){
  if (Stn_name == 'S1'){
      stn_lon  = 145
      stn_lat  = 30
  } else if(Stn_name == 'K2'){
      stn_lon  = 160 
      stn_lat  = 47
  } else if(Stn_name == 'HOT'){
      stn_lon  = -158
      stn_lat  = 22.75
  }
  
  cff  = getdata_station('Chl_C',stn_lon,stn_lat)
  cff$data = cff$data[,-1]
  cff$data = 12/cff$data
  cff$time = 360/length(cff$time)*(cff$time-.5)
  return(cff)
}

pdf('CChl_mu_mod.pdf', width=5,height=6,paper='a4')
   op <- par(font.lab = 1,
                family ="serif",
                mar    = c(4,4,1,3),
                mgp    = c(2.3,1,0),
                mfrow  = c(3,2)) 
j=0
for (Var in c('mu','g','C_Chl')){
  if (Var == 'mu'){
    Ymax = 1.5
  } else if (Var == 'g'){
    Ymax = 1.5
  } else if (Var == 'C_Chl'){
    Ymax = 500
  }
  for (Stn in c('S1','K2')){
    XLab    = ''
    YLab    = ''
    par(mar=c(2,2,2,.2))
    #Write out Ylab:
    if (Var == 'mu' && Stn =='S1'){
       YLab = expression(paste('Growth rate ('*d^-1*')'))
       par(mar=c(2,4,2,.2))
    }else if(Var == 'g'     && Stn=='S1'){
       YLab = expression(paste('Grazing rate ('*d^-1*')'))
       par(mar=c(2,4,2,.2))
    }else if(Var == 'C_Chl' && Stn=='S1'){
       YLab = expression(paste("C:Chl "*' (g '*g^-1*')'))
       par(mar=c(4,4,1,.2))
    }

    if (Var == 'C_Chl') {
       XLab='Date of the year'
    }
    if (Var == 'C_Chl' && Stn=='K2') {
      par(mar=c(4,2,1,.2))
      Ymax=120
    }

  #Obtain MLD data:
    MLD=get_MLD(Stn)
    j  =j+1
    ii =0
    N  =length(MLD)
    #Plot the general frame:
    plot(1:N,1:N,type='n',ylim=c(0,Ymax),
         xlab=XLab,ylab=YLab,xaxt='n')
    axis(1, at=seq(0,360,by=60))  #plot the axis every 2 months
    mtext(letters[j],adj=0,cex=1)
    if (j %in% 1:2) mtext(Stn,adj=.5)
    for (model in c('EFTsimple','Geidersimple')){
       ii      =  ii+1
       filedir =  paste('~/working/FlexEFT1D/DRAM_0.2/',
                          Stn,'/',model,sep='')
       setwd(filedir)
 
       #Get modeled Chl data
       if (Var == 'mu'){

          muvar = paste('muN', formatC(1,width=2,flag= ''),sep='')
           dat  = getData(muvar)
       }else if (Var == 'g'){
          muvar = paste('Gra', formatC(1,width=2,flag= ''),sep='')
           dat  = getData(muvar)
       }else if (Var == 'C_Chl'){
          if (model == 'Geidersimple'){
             CHL  = getData('CHL_T')
             PHY  = paste('PHY', formatC(1,width=2,flag= ''),sep='')
             PHY  = getData(PHY)
             days = PHY$days
            depth = PHY$depth
             CHL  = CHL$data
             PHY  = PHY$data
           #Get Q0N (N:C ratio) from bestpar
             bestpar <- read.table('bestpar',skip=5)
             kk      <- which(bestpar[,1] == 'Q0N')
             Q0N     <- bestpar[kk,4]
               #Calculate Theta:
             dat  = CHL/(PHY/Q0N)
             dat  = list(days=days,depth=depth,data=dat)
          } else{
             THE  = paste('THE', formatC(1,width=2,flag= ''),sep='')
             dat  = getData(THE)
          }
       }
        days = dat$days
       depth = dat$depth
        dat  = dat$data
        d_per_y = 360
        #Get the data of the final year
        w       = (nrow(dat)-d_per_y+1):nrow(dat)
        dat     = dat[w,]

        s_dat   = numeric(length(MLD))
        for (k in 1:length(MLD)){
             w    = which(depth >= -abs(MLD[k]))
             xff  = mean(as.numeric(dat[k,w]))
             if (Var=='C_Chl') xff=12/xff
          s_dat[k]= xff
        }

        setwd('~/Working/FlexEFT1D/')  
        lines(1:length(s_dat),s_dat, lwd=1,col=ii+1)

        #Plot remote sensing estimates of C:Chl ratio on the figure
        if (Var == 'C_Chl'){
           theta = get_theta(Stn)
           points(theta$time,theta$data,pch=2)
        }
        if (Var == 'mu' && Stn == 'S1') {
           legend('topright',c('NPZDFlex','NPZDChl'),col=2:3,lty=1)
        }
    }
  }
}
dev.off()
