#plot histograms of EFTdiscrete model at 4 seasons 
setwd('~/Working/FlexEFT1D/')  
source('~/Working/FlexEFT1D/Rscripts/get_MLD.R')
source('~/Working/FlexEFT1D/Rscripts/getData.R')
source('~/Working/FlexEFT1D/Rscripts/interp.R')
source('~/Working/FlexEFT1D/Rscripts/MLDmean.R')
NPHY    = 20
PMU_min = log(10*pi/6*0.5**3)
PMU_max = log(10*pi/6*60**3)
PMU     = seq(PMU_min,PMU_max,length.out=NPHY)

Vol2D = function(PMU)(exp(PMU)*6/10/pi)**(1/3)
#For QN, Theta, 

for (Var in c('PHY','QN','THE','CHL')){
   fname = paste(Var,'_hist.pdf',sep='')
   pdf(fname, width=5,height=8,paper='a4')
      op <- par(font.lab = 1,
                   family ="serif",
                   mar    = c(2,2,2,0),
                   mgp    = c(2.3,1,0),
                   oma    = c(4,4,0,0),
                   mfcol  = c(4,2)) 
   
   for (Stn in c('S1','K2')){
     #Obtain MLD data:
     MLD  = get_MLD(Stn)
     N    = length(MLD)
     model= 'EFTdiscrete'
     filedir  <- paste('~/working/FlexEFT1D/DRAM_0.3/',
                              Stn,'/',model,sep='')
     setwd(filedir)
     phy = matrix(NA, nr = NPHY, nc = N)
     QN  = phy
     THE = QN
     if (Var == 'QN'){
        ndigit = 3
     } else{
        ndigit = 2
     }
     for (i in 1:NPHY){
       if (Var == 'CHL'){
         Phyvar   = paste('PHY', formatC(i,width=2,flag= ''),sep='')
         phy[i,]  = MLDmean(Phyvar,MLD)
         QNvar    = paste('QN',  formatC(i,width=3,flag= ''),sep='')
         QN[i,]   = MLDmean(QNvar,MLD)
         THEvar   = paste('THE',  formatC(i,width=2,flag= ''),sep='')
         THE[i,]  = MLDmean(THEvar,MLD)
         phy[i,]  = phy[i,]/QN[i,]*THE[i,]
       } else{
         modvar   = paste(Var, formatC(i,width=ndigit,flag= ''),sep='')
         phy[i,]  = MLDmean(modvar,MLD)
       }
     }
     seasons=c('Winter','Spring','Summer','Fall')
     for (i in 1:4){
         #Winter (Jan-Mar average)
         Phy_ = phy[,(1+(i-1)*N/4):(i*N/4)]
         Phy_ = apply(Phy_, 1, mean)
         if (Var == 'THE') Phy_ = 12/Phy_
         #Plot histogram:
         barplot(Phy_,names.arg =round(Vol2D(PMU),1),beside=T)
         mtext(seasons[i],adj=.7)
         if (i == 1) mtext(Stn,adj=0.2)
     }
     setwd('~/Working/FlexEFT1D/')  
   }
   mtext('Phytoplankton size class (µm)',outer=T, adj=.5,side=1)
   if (Var == 'PHY'){
      mtext('Biomass (µM)',outer=T, adj=.5,side=2)
   } else if (Var == 'QN'){
      mtext(expression(paste(Q[N]*' (mol:mol)')),outer=T, adj=.5,side=2)
   } else if (Var == 'THE'){
      mtext('C:Chl (g:g)',outer=T, adj=.5,side=2)
   } else if (Var == 'CHL'){
      mtext('Chl a (µg/L)',outer=T, adj=.5,side=2)
   }
   dev.off()
}
