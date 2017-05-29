#Plot the model output of surface Chl:
source('~/Working/FlexEFT1D/Rscripts/getData.R')
setwd('~/Working/FlexEFT1D/')  
source('Rscripts/interp.R')
source('Rscripts/get_MLD.R')
source('Rscripts/Hz.R')

get_Deu <- function(Stn,model,th=1E-3){
    #Par data: 
    setwd(paste('~/Working/FlexEFT1D/DRAM_0.2/',Stn,'/',model,'/',sep=''))
    dat     = getData('PAR')
    days    = dat$days
      depth = dat$depth
       dat  = dat$data

    d_per_y = 360
    #Get the data of the final year
    w       = (nrow(dat)-d_per_y+1):nrow(dat)
    dat     = dat[w,]
    Deu     = numeric(d_per_y)
    for (i in 1:d_per_y){
        cff = as.numeric(dat[i,])
        cff = cff/cff[length(cff)]
        Deu[i] = depth[which.max(cff >= th) ]
    }
    return(Deu)
}

pdf('ML_mod_obs.pdf', width=5,height=6,paper='a4')
   op <- par(font.lab = 1,
                family ="serif",
                mar    = c(4,4,1,3),
                mgp    = c(2.3,1,0),
                mfrow  = c(3,2)) 
j=0
for (Var in c('TIN','CHL','NPP')){
  for (Stn in c('S1','K2')){

  #Obtain MLD data:
    MLD=get_MLD(Stn)
    j=j+1
    #Read observational data
    dat     = paste('~/Working/FlexEFT1D/',Stn,'/',Stn,'_',Var,'.dat',sep='')
    dat     = read.table(dat,header=T)
    
    #Calculate the average data within the MLD:
    DOY_uni = unique(sort(dat$DOY))
    N       = length(DOY_uni)
    dat2    = numeric(N)

    if (Var == 'NPP') Deu = abs(get_Deu(Stn,'Geidersimple'))
    for (i in 1:N){
       DOY_ = DOY_uni[i]
       MLD_ = MLD[DOY_]
       dat_ = dat[dat$DOY==DOY_,]
       xdep = dat_$Depth
       if (Var == 'NPP'){
       #Calculate integrated NPP (from surface to 0.1% light depth):
          #Get the depth of 0.1% of I0
          xff = 0
         Deu_ = Deu[DOY_]
          Y   = dat_[,3]
          if (max(xdep) < Deu_) {
             xdep = c(xdep,Deu_)
              Y   = c(Y,0)
          }
          for (k in 1:nrow(dat_)){
              xff=xff+(xdep[k+1]-xdep[k])*(Y[k]+Y[k+1])/2
          }
          dat2[i]=xff

       } else{
       #If not NPP:
            dat2_ = dat_[dat_$Depth <= abs(MLD_),]
          if (nrow(dat2_) < 1) {
             dat2[i]= dat_[which.min(xdep),3]
          } else{
             dat2[i]= mean(dat2_[,3])
          }
       }
    }
    ymax    = max(dat2)
    XLab    = ''
    YLab    = ''
    par(mar=c(2,2,2,.2))
    #Write out Ylab:
    if (Var == 'TIN' && Stn =='S1'){
       YLab = expression(paste('TIN (µmol '*L^-1*')'))
       par(mar=c(2,4,2,.2))
    }else if(Var == 'CHL' && Stn=='S1'){
       YLab = expression(paste("Chl "*italic(a)*' (µg '*L^-1*')'))
       par(mar=c(2,4,2,.2))
    } else if(Var == 'NPP' && Stn=='S1'){
       YLab = expression(paste('NPP  (mg C '*m^-2*' '*d^-1*')'))
       par(mar=c(4,4,2,.2))
    }


    if (Var == 'NPP') {
       XLab='Date of the year'
    }
    
    if (Var == 'NPP' && Stn=='K2') par(mar=c(4,2,2,.2))
    #Plot data points:
    plot(DOY_uni,dat2,
      xlab=XLab,
      ylab=YLab,
      xlim=c(0,360),ylim=c(0,ymax),cex=.5,pch=16,xaxt='n')
    axis(1, at=seq(0,360,by=60))
    mtext(letters[j],adj=0,cex=1)
    if (j %in% 1:2) mtext(Stn,adj=.5)

    ii=0
    for (model in c('EFTsimple','Geidersimple')){
       ii=ii+1
       filedir  <- paste('~/working/FlexEFT1D/DRAM_0.2/',
                          Stn,'/',model,sep='')
       setwd(filedir)
 
       #Get modeled Chl data
       if (Var == 'CHL'){
        Chl  = getData('CHL_T')
       }else if (Var == 'TIN'){

        Chl  = getData('NO3')

       }else if(Var == 'NPP'){
        Chl  = getData('NPP_T')
       }
        days = Chl$days
       depth = Chl$depth
        Chl  = Chl$data
        d_per_y = 360
        #Get the data of the final year
        w       = (nrow(Chl)-d_per_y+1):nrow(Chl)
        Chl     = Chl[w,]

        s_Chl   = numeric(length(MLD))

        if (Var == 'NPP'){
        #Calculate integrated values within MLD
           HZ =Hz()
           for (k in 1:length(MLD)){
               s_Chl[k]=sum(HZ*Chl[k,])
           }
        }else{
           
        #calculate the average data within the MLD
           for (k in 1:length(MLD)){
              w    = which(depth >= -abs(MLD[k]))
              xff  = mean(as.numeric(Chl[k,w]))
           s_Chl[k]= xff
           }
        }
        setwd('~/Working/FlexEFT1D/')  
        lines(1:length(s_Chl),s_Chl, lwd=1,col=ii+1)
        if (Var == 'TIN' && Stn == 'S1') legend('topright',c('NPZDFlex','NPZDChl'),col=2:3,lty=1)
    }
  }
}
dev.off()
