#Read data files from 1D model output:
library(plot3D)
plot_1D = function(var){
   file = paste(var,'.out',sep='')
   data = read.table(file,header=F)
   days = as.numeric(as.character(data[2:nrow(data), 2]))
   depth= as.numeric(data[1,3:ncol(data)])
   data = data[2:nrow(data), 3:ncol(data)]
   image2D(as.matrix(data), x=days, y=depth, 
           xlab='Days',ylab='Depth (m)')
}
