getData <- function(var){
   var    <- as.character(var)
   file   <- paste(var,'.out',sep='')
   data   <- read.table(file,header=F)
   days   <- as.numeric(as.character(data[2:nrow(data), 2]))
   depth  <- as.numeric(data[1,3:ncol(data)])
   data   <- data[2:nrow(data), 3:ncol(data)]
   return(list(days=days,depth=depth,data=data))
}

