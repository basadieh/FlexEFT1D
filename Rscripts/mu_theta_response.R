Stn   = 'S1'
model = 'Geidersimple'

#Read best parameters:
filedir = paste('~/working/FlexEFT1D/DRAM_0.2/',
                       Stn,'/',model,sep='')
setwd(filedir)
bestpar = read.table('bestpar',skip=5)

