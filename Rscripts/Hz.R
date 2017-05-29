Hz <- function(hmax=250, thetaS=2, nlev=30){
  Z_w=numeric(nlev+1)
  Z_r=numeric(nlev)
  Hz =numeric(nlev)
  Z_w[1] = -hmax
  for( i   in 1:nlev){
      sc_r    = ((i-nlev) - 0.5)/nlev
     C_sig    = sinh(thetaS*sc_r)/sinh(thetaS)      # -1 < C_sig < 0
     Z_r[i]   = C_sig*hmax
      sc_r    = (i-nlev)/nlev
     C_sig    = sinh(thetaS*sc_r)/sinh(thetaS)      # -1 < C_sig < 0
     Z_w[i+1] = C_sig*hmax
      Hz[i]   = Z_w[i+1] - Z_w[i]
  }
  return(Hz)
}

