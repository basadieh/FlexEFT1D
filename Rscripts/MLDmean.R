#Get the data of the final year and calculate the mean within MLD:
#Chl: a matrix with ncol = depth
source('~/Working/FlexEFT1D/Rscripts/getData.R')
MLDmean <- function(Var,MLD){  
   d_per_y = 360
   Chl     = getData(Var)
   days    = Chl$days
   depth   = Chl$depth
   Chl     = Chl$data
   w       = (nrow(Chl)-d_per_y+1):nrow(Chl)
   Chl     = Chl[w,]
   s_Chl   = numeric(length(MLD))

   #calculate the average data within the MLD
   for (k in 1:length(MLD)){
      w    = which(depth >= -abs(MLD[k]))
      xff  = mean(as.numeric(Chl[k,w]))
   s_Chl[k]= xff
   }
   return(s_Chl)
}

