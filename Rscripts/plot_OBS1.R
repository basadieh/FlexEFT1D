
#Plot observational data:
source('Rscripts/surface_Chl.R')
source('Rscripts/plot_vertical_2D.R')

plot_vertical_2D('TIN','S1')
plot_vertical_2D('Chl','S1')
plot_vertical_2D('NPP','S1')
plot_vertical_2D('TIN','K2')
plot_vertical_2D('Chl','K2')
plot_vertical_2D('NPP','K2')

#Plot 
