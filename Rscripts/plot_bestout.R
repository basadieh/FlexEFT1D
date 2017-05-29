plot_bestout <- function(Stn,model_ID){
     source('~/Working/FlexEFT1D/Rscripts/plot_1D.R')
     #Plot best model outputs:
     filedir  <- paste('~/working/FlexEFT1D/DRAM_0.2/',
                        Stn,'/',model_ID,sep='')
     setwd(filedir)
     mulim  <- c(0,2.5)
     LNOlim <- c(0,1)
     SIlim  <- c(0,1)
     Thelim <- c(0.05,0.6)
     if (Stn == 'K2'){
        Templim <- c(1,12)
         Akslim <- c(0,0.05)
         NO3lim <- c(10,36)
         PHYlim <- c(0.01,0.4)
         ZOOlim <- c(0.01,3.2)
         CHLlim <- c(0.01,1)
         NPPlim <- c(0.1,20)
     } else if(Stn == 'S1'){
        Templim <- c(16,28)
         Akslim <- c(1E-5,0.5)
         NO3lim <- c(0.1,4.5)
         PHYlim <- c(0.01,0.4)
         ZOOlim <- c(0.01,1.0)
         CHLlim <- c(0.01,1)
         NPPlim <- c(0.1,10)
     }
     pdffile <- paste(Stn,model_ID,'_mod_bestout.pdf',sep='')
     pdf(pdffile, width=5,height=4,paper='a4')
        op <- par(font.lab = 1,
                     family ="serif",
                     mar    = c(4,4,1,3),
                     mgp    = c(2,1,0),
                     mfrow  = c(1,1))  
        #Plot temperature:
        plot_1D('Temp',  title = 'Temperature (ºC)',Templim,T)

        #Plot vertical diffusivity:
        plot_1D('Aks',  
        title = expression(paste('Vertical diffusivity ('
                                *m^2*' '*s^-1*')',sep='')),Akslim,T)
        #Plot NO3:
        plot_1D('NO3',   title = 'NO3 (µM)', NO3lim,T)
      
        #Plot PHY:
        plot_1D('PHY 1', title = 'PHY (µM)',PHYlim, T)
        #Plot ZOO:
        plot_1D('ZOO',title='ZOO (µM)',ZOOlim, T)
        #Plot CHL:
        plot_1D('CHL_T', 
        title = expression(paste('Chl (µg '*L^-1*' )')),CHLlim, T)

        #Plot CHL:C ratio:
        plot_1D('The 1', 
        title = expression(paste('Chl:C (gChl '*molC^-1*' )')),Thelim ,T)

        plot_1D('NPP_T', 
        title = expression(paste('NPP (µgC '*L^-1*' '*d^-1*' )')),NPPlim,T)

        muvar    <- paste('muN', formatC(1,width=2,flag= ''),sep='')
        plot_1D(muvar,   
        title = expression(paste('Net growth rate ('*d^-1*' )')),mulim, T)

        plot_1D('LNO 1',title='Nutrient limitation index',LNOlim,T)
        plot_1D('SI  1',title='Light limitation index',SIlim,T)
     dev.off()
}


