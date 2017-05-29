#Compare PHY, ZOO, and DET between the two models within ML
setwd('~/Working/FlexEFT1D/')  
source('Rscripts/get_MLD.R')
source('~/Working/FlexEFT1D/Rscripts/interp.R')
source('~/Working/FlexEFT1D/Rscripts/MLDmean.R')
pdf('PHY_ZOO.pdf', width=5,height=7,paper='a4')
 op <- par(font.lab = 1,
             family ="serif",
             mar    = c(2,2,1.5,0),
             mgp    = c(2.3,1,0),
             mfrow  = c(3,2),
             oma    = c(4,4,0,0)) 
j=0
for (Var in c('PHY 1','ZOO','DET')){
  if (Var == 'PHY 1'){
    Ymax = .5
  } else if (Var == 'ZOO'){
    Ymax = 1.2
  } else if (Var == 'DET'){
    Ymax = 1.5
  }
  for (Stn in c('S1','K2')){
    XLab = ''
    YLab = ''

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
 
       s_dat = MLDmean(Var, MLD)
       setwd('~/Working/FlexEFT1D/')  
       lines(1:length(s_dat),s_dat, lwd=1,col=ii+1)

       #Plot legend:
       if (Var == 'PHY 1' && Stn == 'S1') {
          legend('topright',c('NPZDFlex','NPZDChl'),col=2:3,lty=1)
       }
    }
  }
}
mtext('Date of the year',side=1,adj=.5,outer=TRUE,line=.5)
vars=c('PHY','ZOO','DET')
dis =c(2.7,1.5,.3)
for (i in 1:length(vars)){
    txt  <- bquote(.(vars[i]) ~ ' (mmol '*m^-3*')')
    mtext(txt,side=2,adj=dis[i]/3,outer=TRUE)
}
par(op)
dev.off()
