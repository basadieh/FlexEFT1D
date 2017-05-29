#Test the subroutine of calculating light attenuation:
vertical_light <- function(I_0,nlev,z_w,Chl){
          dyn.load('FlexEFT.so')
          PAR <- .Fortran('Calculate_PAR',
                      I_0  = as.numeric(I_0),
                      nlev = as.integer(nlev),
                      z_w  = as.numeric(z_w),
                      Chl  = as.numeric(Chl),
                      PAR  = numeric(nlev))$PAR
          return(PAR)}

vertical_light(I_0 = 200, nlev = 5, z_w = rev(c(0, -2, -4, -6, -8, -10)), 
               Chl = c(0.1,0.3,1,.2,.1))


#Test setup_grid

thetaS = 2
nlev   = 10
hmax   = 200
setup_grid <- function(){
          dyn.load('FlexEFT.so')
          PAR <- .Fortran('setup_grid',
                      thetaS  = as.numeric(thetaS),
                      nlev    = as.integer(nlev),
                      hmax    = as.numeric(hmax),
                      Z_r     = numeric(nlev),
                      Z_w     = numeric(nlev+1),
                      Hz      = numeric(nlev))$PAR

          return(list(Z_r = Z_r, Z_w = Z_w, Hz = Hz))
          }
