get_MLD <- function(Stn_name,var='Aks',useTaketo=FALSE){
    if (Stn_name == 'K2') useTaketo=TRUE
    outfile  =  paste('~/Working/FlexEFT1D/',
                      Stn_name,'/',Stn_name,'_',var,'.dat',sep='')
    #Write out timefile:
    timefile =  paste('~/Working/FlexEFT1D/',
                      Stn_name,'/',Stn_name,'_',var,'_time.dat',sep='')

    pdffile  = paste('~/Working/FlexEFT1D/',
                     Stn_name,'/',Stn_name,'_',var,'.pdf',sep='')

    if (useTaketo && var == 'Aks'){
     timefile= paste('~/Working/FlexEFT1D/',Stn_name,'/',
                          Stn_name,'_Aks_timeT.dat',sep='')

     outfile = paste('~/Working/FlexEFT1D/',Stn_name,'/',Stn_name,'_AksT.dat',sep='')

    }
    time = read.table(timefile,header=T)
    data = read.table(outfile, header=T)

    time  = as.numeric(time)
    M     = length(time)
    depth = data[,1]
    data1 = as.matrix(t(data[,2:ncol(data)]))

    if (var == 'Aks'){
       #Calculate MLD based on threshold of Kv (1E-4 m2/s):
       MLD <- numeric(M)
       for (i in 1:M){
           MLD[i] <- interp(x=depth,y=data[,i+1],Xth=1E-4)
       } 
    }
    return(MLD)
}

