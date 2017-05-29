setwd('~/Working/FlexEFT1D/')  
source('~/Working/FlexEFT1D/Rscripts/getData.R')

for (Var in c('TIN','CHL','NPP','muN 1','C_Chl','PHY 1','ZOO','DET')){
  fname = paste(Var,'_Vertical_mod_obs.pdf',sep='')
  pdf(fname, width=5,height=6,paper='a4')
  op <- par(font.lab = 1,
               family ="serif",
               mar    = c(2,2,1.5,0),
               mgp    = c(2.3,1,0),
               oma    = c(4,4,0,0),
               mfcol  = c(4,2)) 
  j=0
  for (Stn in c('S1','K2')){
    #Read observational data
    if (Var %in% c('TIN','CHL','NPP')){
       dat = paste('~/Working/FlexEFT1D/',Stn,'/',Stn,'_',Var,'.dat',sep='')
       dat = read.table(dat,header=T)
    }
    #Average into 4 seasons
    DOYs = seq(0,360,length.out=5)
    seasons <- c('Winter','Spring','Summer','Fall')
    for (i in 1:4){
      j    = j+1
      if (Var %in% c('TIN','CHL','NPP')){
        dat_ = dat[dat$DOY > DOYs[i] & dat$DOY <= DOYs[i+1],]
        plot(dat_[,3],-dat_$Depth,xlab='',ylab='',pch=16,cex=.5)
      } else{
        Depth = -250:-1
        if (Var == 'muN 1'){
           Ymax  = 1.5
        } else if(Var == 'C_Chl'){
           if (Stn == 'S1'){
              Ymax  = 500
           } else if(Stn == 'K2'){
              Ymax  = 120
           }
        } else if(Var == 'PHY 1'){
           Ymax=.5
        } else if(Var == 'ZOO'){
           Ymax=1.2
        } else if(Var == 'DET'){
           Ymax=1.5
        }
        plot(rep(Ymax,length(Depth)),Depth,
             xlim=c(0,Ymax),
             xlab='',ylab='',type='n')
      }
      mtext(paste(letters[j],')',seasons[i]),adj=0)
      if (i==1) mtext(Stn,adj=.7)
      ii=1
      for (model in c('EFTsimple','Geidersimple')){
         ii=ii+1
         filedir  <- paste('~/working/FlexEFT1D/DRAM_0.2/',
                            Stn,'/',model,sep='')
         setwd(filedir)
         #Get modeled data
         if (Var == 'CHL'){
          Chl  = getData('CHL_T')
         }else if (Var == 'TIN'){
          Chl  = getData('NO3')
         }else if(Var == 'NPP'){
          Chl  = getData('NPP_T')
         }else if(Var == 'C_Chl'){
            if (model == 'EFTsimple'){
               Chl  = getData('THE 1')
            }else{
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
             dat     <- CHL/(PHY/Q0N)
             Chl     <- list(days=days,depth=depth,data=dat)
            }
         }else{
          Chl  = getData(Var)
         }
          days    = Chl$days
         depth    = Chl$depth
          Chl     = Chl$data
          d_per_y = 360
          #Get the data of the final year
          w       = (nrow(Chl)-d_per_y+1):nrow(Chl)
          Chl     = Chl[w,]

          #Read model results:
          mdoy  = 1:nrow(Chl)
          k     = which(mdoy > DOYs[i] & mdoy <= DOYs[i+1])
          Chl_  = Chl[k,]
          #Calculate quantiles:
          cff   = sapply(1:ncol(Chl_),function(k)quantile(Chl_[,k], probs=c(0.025,0.5,0.975)))
          setwd('~/Working/FlexEFT1D/')  
          if (Var == 'C_Chl') cff = 12/cff
          matlines(t(cff), depth, lty=c(2,1,2), lwd=c(1,1.5,1), col=ii)
      }
      if (j==1) {
         if (Var == 'TIN' | Var == 'C_Chl'){
            legend('topright',c('NPZDFlex','NPZDChl'),col=2:3,lty=1)
         }else{
            legend('bottomright',c('NPZDFlex','NPZDChl'),col=2:3,lty=1)
         }
      }
    }
  }   

  if (Var == 'TIN'){
    Varname  <- bquote(.(Var)~' (mmol '*m^-3*')')
  }else if (Var == 'CHL'){
    Varname  <-  bquote(.(Var)~' (mg '*m^-3*')')
  }else if (Var == 'NPP'){
    Varname  <- bquote(.(Var)~' (mg C '*m^-3*' '*d^-1*')')
  }else if (Var == 'muN 1'){
    Varname  <- bquote('µ ( '*d^-1*')')
  }else if (Var == 'C_Chl'){
    Varname  <- expression(paste("C:Chl "*' (g '*g^-1*')'))
  }else if (Var == 'PHY 1'){
    Varname  <- expression(paste("PHY "*'(mmol '*m^-3*')'))
  }else{
    Varname  <- bquote(.(Var) ~ ' (mmol '*m^-3*')')
  }

  mtext(Varname,    side = 1, outer=TRUE,line=1)
  mtext('Depth (m)',side = 2, outer=TRUE,line=1)
  dev.off()
}
