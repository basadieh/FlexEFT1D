#Plot mean size and variance for the two size models in MLD
setwd('~/Working/FlexEFT1D/')  
source('~/Working/FlexEFT1D/Rscripts/get_MLD.R')
source('~/Working/FlexEFT1D/Rscripts/getData.R')
source('~/Working/FlexEFT1D/Rscripts/interp.R')
source('~/Working/FlexEFT1D/Rscripts/MLDmean.R')
pdf('Meansize_MLD.pdf', width=5,height=3,paper='a4')
   op <- par(font.lab = 1,
                family ="serif",
                mar    = c(4,4,1,.2),
                mgp    = c(2.1,1,0),
                mfrow  = c(1,2), 
                oma    = c(4,4,0,0)) 

j = 0
for (Stn in c('S1','K2')){
  #Obtain MLD data:
  MLD  = get_MLD(Stn)
  N    = length(MLD)
  Ymax = 60

  par(mar=c(2,2,2,.2))
  XLab=''
  YLab=''

   #Plot the general frame:
    j = j + 1
    plot(1:N,1:N,type='n',ylim=c(0,Ymax),
         xlab=XLab,ylab=YLab,xaxt='n')
    axis(1, at=seq(0,360,by=60))  #plot the axis every 2 months
    mtext(letters[j],adj=0,cex=1)
    ii = 1
    for (model in c('EFTdiscrete', 'EFTcont')){
       ii = ii + 1 #For label color
       filedir  <- paste('~/working/FlexEFT1D/DRAM_0.2/',
                             Stn,'/',model,sep='')
       setwd(filedir)
       PMUavg = MLDmean('PMU',MLD)
       VARavg = MLDmean('VAR',MLD)
       #Calculate the mean and 95% CI of size
       PMU = matrix(NA, nr = length(PMUavg), nc = 3)
       PMU[,1] = PMUavg
       PMU[,2] = PMUavg - 2 * sqrt(VARavg)
       PMU[,3] = PMUavg + 2 * sqrt(VARavg)
       PMU     = (exp(PMU)*6/10/pi)**(1/3)
       setwd('~/Working/FlexEFT1D/')  
       matlines(1:nrow(PMU),PMU, lwd=c(1,.5,.5),lty=c(1,2,2),col=ii)
       if (Stn == 'S1') {
          legend('topleft',c('Discrete','Continuous'),col=2:3,lty=1,cex=.8)
       }
     }
}

mtext('Date of the year',side=1,adj=.5,outer=T)
mtext('Mean size (µm)',  side=2,adj=.5,outer=T)
dev.off()
