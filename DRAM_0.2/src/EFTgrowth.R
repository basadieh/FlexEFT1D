source('~/Working/Adap_uptake/lambert.R')
TEMPBOL <- function(t, Ep=0.4, Tr = 15){
   kb <- 8.62E-5
   t  <- t+273.15
   Tr <- Tr+273.15
   return(exp(-(Ep/kb)*(1/t-1/Tr)))
}

mu <- function(mu0=5,V0=5,PAR=300,NO3=10,Temp=15,aI0=0.15,A0N,Q0N){
   tf_p <- TEMPBOL(Temp)
   Qs   <- Q0N/2
 mu0hat <- tf_p*mu0
  V0hat <- tf_p*V0
  # Initial slope of P-I curve
  aI    <- aI0
  # Cost of photosynthesis
  RMchl0 <- 0.1
  RMchl <- tf_p*RMchl0
  # Threshold irradiance and RMchl is set temperature dependent
  zetaChl <- 0.8
  I_zero <- zetaChl*RMchl/aI  
  
  A0hat <- tf_p*A0N
    fA  <- 1/(1 + sqrt(A0hat*NO3/V0hat) ) 
  VNhat <- (1-fA)*V0hat*fA*A0hat*NO3/((1-fA)*V0hat+fA*A0hat*NO3) 

  if( PAR > I_zero ) {
    larg1 = exp(1 + aI*PAR/(mu0hat*zetaChl))
    larg  = (1 + RMchl/mu0hat)*larg1   
    
     W_Y  = WAPR(larg,0,0)

   ThetaHat = 1/zetaChl + (1- W_Y)*mu0hat/(aI * PAR)

    SI = 1 - max(exp(-aI*PAR*ThetaHat/mu0hat),0)
    mu0hatSI = mu0hat*SI
   
    muIhat   = mu0hatSI-(mu0hatSI+RMchl)*zetaChl*ThetaHat
   
    zetaN  = 0.6
    ZINT   = Qs*(muIhat/VNhat + zetaN)
           
    fV = (-1 + sqrt(1 + 1/ZINT))*Qs*muIhat/VNhat
    fV = max(fV,0.01)
    
    }else{
    # Under the conditions of no light:
       ThetaHat      = 0.01 
       ZINT          = Qs*zetaN
       fV            = 0.01
       muIhat        = 0
    }
     # Nutrient limitation index:
     Lno3 =1/(1 + sqrt(1 +1/ZINT)) 
     QN   = Qs/Lno3
    
     if (PAR > I_zero) {
   #Net growth rate (d-1) of phytoplankton at the average size
       muNet = muIhat*(1-2*Lno3)
     }else{
       muNet = 0
     }
     return(muNet)

}

A0N <- seq(0.01,1,0.01)
mu0 <- c(1,5,10,20)
muNet <- matrix(NA, nr = length(A0N), nc=length(mu0))
for (i in 1:length(A0N)){
    for (j in 1:length(mu0)){
         muNet[i,j] <- mu(A0N = A0N[i],mu0=mu0[j],Q0N=0.01)
    }
}
plot(A0N,muNet[,length(mu0)],ylim=c(0,8),type = 'n')
for (i in 1:length(mu0)){
    points(A0N,muNet[,i],col=i,type='l')
}
