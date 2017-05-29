PROGRAM main
  implicit none
  character(LEN=2 ), parameter :: Stn = 'K2'
  ! Input files
  character(LEN=20), parameter :: Tempfile = trim(Stn)//'_Temp.dat', &
                                   PARfile = trim(Stn)//'_par.dat' , &
                                   Aksfile = trim(Stn)//'_Aks.dat' , &
                                   NO3file = trim(Stn)//'_NO3.dat'  
                                   
  ! !INPUT PARAMETERS for diffusion:
  integer,parameter   :: Dirichlet      = 0
  integer,parameter   :: Neumann        = 1
                  
 
  integer, parameter :: nlev = 40,   &   ! Number of vertical layers
                        NVAR = 7,    &   ! Number of total state variables
                        ! Indices for state variables
                        iNO3 = 1, iPHY = 2, iZOO = 3, iDET = 4, &
                        iPMU = 5, iVAR = 6, iCHL = 7,           &
                        NVsinkterms = 2,                        &
                        Windex(NVsinkterms) = (/iPHY, iDET/),   &
                        ! Number of vertical points of NO3 profile
                        N_NO3  = 14, &  
                        N_Temp = 24, &
                        N_Aks  = 40, &
                        Nout   = 28, &     ! Number of output variables
                        oTemp  = 1,  &
                        oPAR   = 2,  &
                        oAks   = 3,  &
                        oNO3   = 1,  &
                        oPHY   = 2,  &
                        oZOO   = 3,  &
                        oDET   = 4,  &
                        oPMU   = 5,  &
                        oVAR   = 6,  &
                        oCHL   = 7, &
                        oFER   = 8, &
                        oTheta = 9, &
                        oQN    = 10, &
                        omuNet = 11, &
                        oGraz  = 12, &
                        oZ2N   = 13, &
                        oD2N   = 14, &
                        odmudl = 15, &
                        od2mu  = 16, &
                        od2gdl = 17, &
                        ow_p   = 18, &
                  ! Diffusion fluxes
                        oD_NO3 = 19, &
                        oD_PHY = 20, &
                        oD_ZOO = 21, &
                        oD_DET = 22, &
                        oD_PMU = 23, &
                        oD_VAR = 24, &
                        oD_CHL = 25
                        

  ! Output indices:
  character(LEN=5), parameter :: Labelout(Nout)  &
  = (/'Temp ','PAR  ','Aks  ','NO3  ','PHY  ','ZOO  ','DET  ','PMU  ',&
      'VAR  ','CHL  ','FER  ','Theta','QN   ','muNet','Graz ','Z2N  ',&
      'D2N  ','dmudl','d2mu ','d2gdl','w_p  ',                        &
      'D_NO3','D_PHY','D_ZOO','D_DET','D_PMU','D_VAR','D_CHL'/)  

  real(8), parameter :: hmax = 250., &   ! Total water depth
                       thetaS= 2.  , &   ! surface stretching parameter
                       dtsec = 60. , &   ! time step in seconds
                     d_per_s = 86400.,&    ! how many seconds in one day
                       dtdays= dtsec/d_per_s, & !the fraction of a time step in one day
                       wDET  = 2d0         ! Sinking rate of Detritus (m/d)

  integer, parameter :: y_per_s = INT(d_per_s*360), &!how many seconds of one year
                        nsave= INT(d_per_s)/INT(dtsec), &! Timesteps to save
                        maxD = 360*5,  &  ! Number of days to run
                        maxN = maxD*INT(d_per_s)/INT(dtsec)  ! Number of time steps

  real(8),parameter   :: DefaultRelaxTau(nlev)= 1.d15, cnpar = 0.6

  real(8) :: Z_r(1:nlev), Z_w(0:nlev), Hz(nlev),  &  ! Grid variables
             Vars(NVar,nlev),                     &  ! State variables
             NewVars(NVar,nlev) ,                 &  ! New calculated state variables
             obs_NO3(N_NO3  , 2),                 &
             obs_PAR(1      ,13),                 &
             obs_Temp(N_Temp,13),                 &
             obs_Aks(N_Aks  ,13),                 &
             obs_time(12), &                ! observation time indices 
             time, I_0, PAR(nlev), Temp(nlev), Aks(0:nlev),     &
             NO3(nlev,1),                                       &
             VTemp(nlev,12), VAks(0:nlev,12), Varout(100,nlev), &
             ! Sinking rate
             ww(0:nlev,NVsinkterms), Lsour(nlev), Qsour(nlev)             

  integer :: i, it, k, j, time_day

  character(LEN=20)  :: outfile

  ! Setup grid
  call setup_grid(thetaS, nlev, hmax, Z_r, Z_w, Hz)

  ! Initialize state variables
  call Readcsv(NO3file, N_NO3,     2, obs_NO3) 

  ! Interpolate initial NO3:
  call gridinterpol(N_NO3,1,obs_NO3(:,1),obs_NO3(:,2),   &
                    nlev, Z_r, NO3(:,1)) 
  ! Initialize initial NO3:
  Vars(iNO3,:) = NO3(:,1)


  ! Initialize other variables:
  do i = 2,NVAR
     do k = 1,nlev
        Vars(i, k) = 0.1
     enddo
  enddo
  ! Initialize output data:
  Varout(:,:) = 0.0
  do i = 1,NVAR
     do k = 1,nlev
        Varout(i,k)= Vars(i,k)
     enddo
  enddo

  ! Give obs. time indices in seconds:
  do i = 1, 12
     obs_time(i) = (float(i-1)+0.5)*30.*d_per_s
  enddo

  ! Read external PAR data:
  call Readcsv(PARfile, 1,      13,obs_PAR) 

  ! Read external Temp data:
  call Readcsv(Tempfile,N_Temp, 13,obs_Temp) 

  ! Interpolate external Temp data:
  call gridinterpol(N_Temp,12,obs_Temp(:,1),obs_Temp(:,2:13),   &
                       nlev, Z_r, VTemp(:,:)) 

  ! Interpolate external Aks data:
  call Readcsv(Aksfile, N_Aks ,   13,obs_Aks) 

  ! Interpolate external Aks data:
  call gridinterpol(N_Aks,12,obs_Aks(:,1), obs_Aks(:,2:13),   &
                       nlev+1, Z_w, VAks(:,:)) 

  ! Create output files:
  do i = 1, Nout
     outfile = trim(Labelout(i))//'.out'

     ! Create data out file
     open (9+i, file = outfile, status = 'replace')
   
     if (i .eq. oAks) then
        write(9+i, 3) 'Timestep','Days',Z_w
     else
        write(9+i, 1) 'Timestep','Days',Z_r
     endif
  enddo

1 format(A10,1x,A7, <nlev  >(2x,F12.5))
2 format(I10,1x,I7, <nlev  >(2x,F12.5))
3 format(A10,1x,A7, <nlev+1>(2x,F12.5))
4 format(I10,1x,I7, <nlev+1>(2x,F12.5))

  ! For each time step, read in external environmental data
  write(6,*) 'START TIME STEPPING'
  DO it = 1, maxN+1

  ! Calculate current timing (zero is starting time):
     time = float(it-1)*dtsec

  ! Interpolate Temp data throughout the water column:
     call time_interp(int(time), 12, nlev, obs_time, VTemp, Temp)

  ! Interpolate temporal PAR data: 
     call time_interp(int(time), 12,1, obs_time,obs_PAR, I_0)

  ! Convert the unit of par to W m-2:
     I_0 = I_0/0.4

  ! Calculate light field:
     call Calculate_PAR(I_0, nlev, Hz, Vars(iCHL,:), PAR)

  ! Interpolate Aks data throughout the water column:
     call time_interp(int(time), 12, nlev+1, obs_time, VAks, Aks)


  ! Save data to output files:
  if (mod(it, nsave) .eq. 1) then

  ! Calculate model time in days:
     time_day=int(time/d_per_s)

  ! Save Temp data into the output file:
     write(9+oTemp, 2) it, time_day, Temp

  ! Save PAR data into the output file:
     write(9+oPAR, 2) it, time_day, PAR

  ! Save Aks data into the output file:
     write(9+oAks, 4) it, time_day, Aks

  ! Save state variables and diagnostics:
     do i=4,Nout
        write(9+i, 2) it, time_day, Varout(i-3,:)
     enddo

  endif  !==> End of saving results

  ! Biological rhs:
  do k = 1,nlev
     call FlexEFT(Temp(k),PAR(k),dtdays,NVAR,Vars(:,k),Varout(:,k))

     ! Pass the new state variables to Vars
     do j = 1,NVAR
        Vars(j,k)=Varout(j,k)
     enddo
  enddo

  ! Diffusion:
  Lsour(:) = 0.
  Qsour(:) = 0.
  do j = 1,NVAR

     call diff_center(nlev,dtsec,cnpar,1,Hz,Neumann,Neumann,&
        0.,0.,Aks,Lsour,Qsour,DefaultRelaxTau,Vars(j,:),Vars(j,:),NewVars(j,:))

     ! Save diffusion fluxes (normalized to per day)
     Varout(18+j,:) = (NewVars(j,:) - Vars(j,:))/dtdays
 
     ! Update the state variables:
     Vars(j,:)=NewVars(j,:)
  enddo

  ! Sinking:
  ! Initialize sinking rate:
  ww(0,   :) = 0.0
  ww(nlev,:) = 0.0

  do k = 1,nlev-1
     ww(k,1) = Varout(ow_p,k)         !Phytoplankton sinking rate  
     ww(k,2) = -wDET                  !Detritus      sinking rate
  enddo

  do j = 1,NVsinkterms
     call adv_center(nlev,dtdays,Hz,Hz,ww(:,j),  &
                     1,1,0.,0.,6,1,Vars(Windex(j),:))
  enddo

  ! Check whether the values are valid:
  do j = 1,NVAR
     do k=1,nlev
        if( (Vars(j,k) .ne. Vars(j,k)) .OR. (Vars(j,k) .le. 0.) ) then
            write(6,*) 'At time step ',it
            write(6,*) 'The variable ',trim(Labelout(j+3)), 'is invalid at depth ',Z_r(k)
            stop
        endif
     enddo 
  enddo

  ENDDO    ! ==> End of time stepping
  ! Close files:
  do i = 1, Nout
     close (unit=9+i)
  enddo

END PROGRAM main
!========================================================
! !ROUTINE: Interpolate from time-series observation to model time
subroutine time_interp(time, N, nrow, obs_time, obs_data, mod_data)
  ! time: model time
  !    N: number of time-series observations in observational data
  ! nrow: number of vertical points in observation data
  ! obs_time: time of observational data, should be a vector of length N
  ! obs_data: observational data, should be a vector of length N
  ! mod_data: model data as output
  implicit none
  integer, intent(in)    :: N, time, nrow
  real(8), intent(in)    :: obs_time(N), obs_data(nrow,N)
  real(8), intent(out)   :: mod_data(nrow)

  integer, parameter     :: y_per_s = 31104000 
  integer                :: i,j
  real(8)                :: timeMOD, rat

! Get time of the year (mod):
   timeMOD = float(mod(time,y_per_s))

! Deal with the case between two years:
   IF ((timeMOD .lt. obs_time(1))) THEN 

      rat = (timeMOD    -obs_time(N)+float(y_per_s))   &
          / (obs_time(1)-obs_time(N)+float(y_per_s))

      do i=1,nrow
         mod_data(i)=obs_data(i,N)*(1d0-rat)+obs_data(i,1)*rat
      enddo

   ELSEIF (timeMOD .ge. obs_time(N)) THEN

      rat = (timeMOD    -obs_time(N))   &
          / (obs_time(1)-obs_time(N)+float(y_per_s))

      do i=1,nrow
         mod_data(i)=obs_data(i,N)*(1d0-rat)+obs_data(i,1)*rat
      enddo

   ELSE

      do j = 1, N-1
         if ((timeMOD.lt.obs_time(j+1)) .and. (timeMOD.GE.obs_time(j))) then
            
            rat=(timeMOD-obs_time(j))/(obs_time(j+1)-obs_time(j))
       
            do i=1,nrow
               mod_data(i)=obs_data(i,j)*(1d0-rat)+obs_data(i,j+1)*rat 
            enddo
            exit
         endif
      enddo 
 
   ENDIF
   return
end subroutine
!========================================================
! !ROUTINE: Interpolate from observation space to model grid
!
! !INTERFACE:
subroutine gridinterpol(N,cols,obs_z,obs_prof,nlev,model_z,model_prof)
!
! !DESCRIPTION:
!
!  This is a utility subroutine in which observational data, which might
!  be given on an arbitrary, but structured grid, are linearly interpolated and
!  extrapolated to the actual (moving) model grid.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   !N    : number of vertical layers in observation data
   !cols : number of profiles in obs. data
   integer,  intent(in)                :: N,cols
   REAL(8),  intent(in)                :: obs_z(N),obs_prof(N,cols)
   integer,  intent(in)                :: nlev
   REAL(8),  intent(in)                :: model_z(nlev)
!
! !OUTPUT PARAMETERS:
   REAL(8),  intent(out)               :: model_prof(nlev,cols)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i,j,ii
   REAL(8)                   :: rat
!-----------------------------------------------------------------------
!BOC
!  Set surface values to uppermost input value
   do i=nlev,1,-1
      if(model_z(i) .ge. obs_z(N)) then
         do j=1,cols
            model_prof(i,j) = obs_prof(N,j)
         end do
      end if
   end do

!  Set bottom values to lowest input value
   do i=1,nlev
      if(model_z(i) .le. obs_z(1)) then
         do j=1,cols
            model_prof(i,j) = obs_prof(1,j)
         end do
      end if
   end do

!  Interpolate inner values linearly
   do i=1,nlev
         !write(*,*) 'i = ',i
      if ((model_z(i) .lt. obs_z(N)) .and. (model_z(i) .gt. obs_z(1))) then
         ii=0
224      ii=ii+1
         if (obs_z(ii) .le. model_z(i)) goto 224
         ! Debug:
         !write(*,*) 'ii = ', ii
         !write(*,*) 'model_z = ',model_z(i)
         !write(*,*) '  obs_z = ',  obs_z(ii)
         !write(*,*) '  obs_z(ii-1) = ',  obs_z(ii-1)

         rat=(model_z(i)-obs_z(ii-1))/(obs_z(ii)-obs_z(ii-1))

         do j=1,cols
            model_prof(i,j)=(1-rat)*obs_prof(ii-1,j)+rat*obs_prof(ii,j)
         !write(*,*) 'obs_prof(ii-1) = ',obs_prof(ii-1,j)
         !write(*,*) 'obs_prof(ii)   = ',obs_prof(ii  ,j)
         !write(*,*) 'rat            = ',rat
         !write(*,*) 'model_prof(i)  = ',model_prof(i,j)

         end do
      end if
   end do

   return
end subroutine
!EOC

!========================================================
subroutine setup_grid(thetaS, nlev, hmax, Z_r, Z_w, Hz)
  implicit none
  integer, intent(in)    :: nlev
  real(8), intent(in)    :: thetaS, hmax
  real(8), intent(out)   :: Z_r(1:nlev), Z_w(0:nlev), Hz(1:nlev)

  ! Local scratch variable:
  integer                :: i
  real(8)                :: sc_r, C_sig
  
  Z_w(0) = -hmax

  !Following Song and Haidvogel (1994). sinh is the hyperbolic sin function
  do i   = 1,nlev
      sc_r  = (float(i-nlev) - 0.5)/float(nlev)
     C_sig  = sinh(thetaS*sc_r)/sinh(thetaS)      ! -1 < C_sig < 0
     Z_r(i) = C_sig*hmax
      sc_r  = (float(i-nlev))/float(nlev)
     C_sig  = sinh(thetaS*sc_r)/sinh(thetaS)      ! -1 < C_sig < 0
     Z_w(i) = C_sig*hmax
      Hz(i) = Z_w(i) - Z_w(i-1)
  enddo
  return
end subroutine

!========================================================
SUBROUTINE FlexEFT(tC,par,dtdays,NVAR,Varin,Varout)
  implicit none
!INPUT PARAMETERS:
  real(8), intent(in)    :: tC,par,dtdays
  integer, intent(in)    :: NVAR
  real(8), intent(in)    :: Varin(NVAR)
  integer, parameter     :: Nout = 100
  real(8), intent(out)   :: Varout(Nout)
!LOCAL VARIABLES of phytoplankton:
  real(8) :: NO3,PHY,ZOO,DET,PMUPHY,VARPHY,NO31,PHY1,ZOO1,DET1,      &
              PMUPHY1,VARPHY1,CHL,CHL1  !State variables
  real    :: tf_p,tf_z,I_zero,larg   !Environmental variables
  real(8) :: PMU,VAR,PMU1,VAR1,X,B,dXdl,dBdl,KFe,Fe,dmu0hatdl, d2mu0hatdl2
  real(8) :: V0hat,Kn,Lmin,A0N,A0hat,dA0hatdl,d2A0hatdl2,fA
  real(8) :: VNhat,dVNhatdl,d2VNhatdl2, fN, d2fNdl2 ! Nutrient uptake variables
  real(8) :: mu0hat,muIhat,mu0hatSI,dmu0hatSIdl,d2mu0hatSIdl2
  real(8) :: dmuIhatdl,d2muIhatdl2! Growth rate variables
  real(8) :: fV,dfVdl,d2fVdl2  !fV
  real(8) :: ZINT,dZINdl,d2ZINdl2 !ZINT
  real(8) :: aI,SI,dSIdl,d2SIdl2         !Light dependent component
  real(8) :: RMchl,Theta,ThetaHat,dThetaHatdl,d2ThetaHatdl2 !Chl related variables
  real(8) :: QN,Qs,dQNdl,d2QNdl2  ! cell quota related variables
  real(8) :: larg1,w_p,dmu0hat_aIdl,d2mu0hat_aIdl2,dlargdl,                       &
          d2largdl2,W_Y,dWYYdl,daI_mu0hatdl,d2aI_mu0hatdl2,                       &
          d2wpdl2  
  real(8) :: dV0hatdl,d2V0hatdl2,muNet,dmuNetdl,d2muNetdl2,ThetaAvg,QNavg,dwdl,   &
              d2wdl2,w_pAvg,dgdlbar,d2gdl2bar
  real(8),parameter :: thetamin = 0.01, DETmin = 1D-2,DETmax = 1D5                &
          ,QNmin = 0.01,fVmin   = 0.01, Lref   = 0.5,    eps = 1D-4               &
          ,PMUmax= 20., VARmax  = 100., Penfac = 1.0,    pi  = 3.14159265358979   &
          ,Ep    = 0.4, alphaQ  = -0.18,alphaI = -0.13,  alphamu = 0.             &
          ,betamu= 0. , Femin   = 0.02, mu0    = 5.0,    aI0 = 0.25               &
          ,A0N0  = 5.0,  alphaA = -0.3, V0N    = 5.0, alphaV = 0.2                &
          ,K0N   = 0.17, alphaK = 0.27, K0Fe   = 0.8,alphaFe = 0.14               &
          ,alphaG = 1.1, zetaChl= 0.8 , zetaN  = 0.6,RMchl0  = 0.1                &
          ,mp     = 0.02,Q0N    = 0.1,  w_p0   = -0.05,alphaW = 0.0,              &
           w_d    = -5., Wmax   = 50.0, gmax   = 1.0,  kp     = 0.5,              &
           rdn    = 0.2, Ez     = 0.6,  mz     = 0.05
  integer, parameter :: nutrient_uptake = 2, grazing_formulation = 3   
  real(8), external :: temp, WAPR, ScaleTrait, PenAff, grazing
  !Declarations of zooplankton:
  real(8) :: Cf,aG,INGES,RES,EGES,gbar
  real(8),parameter :: unass = 0.24, GGE = 0.3333
  integer :: i,j,gform
  integer, parameter :: iNO3=1,iPHY=2,iZOO=3,iDET=4,iPMU=5,iVAR=6,iCHL=7,&
                        iFER=8,      &
                        oNO3   = 1,  &
                        oPHY   = 2,  &
                        oZOO   = 3,  &
                        oDET   = 4,  &
                        oPMU   = 5,  &
                        oVAR   = 6,  &
                        oCHL   = 7,  &
                        oFER   = 8,  &
                        oTheta = 9,  &
                        oQN    = 10, &
                        omuNet = 11, &
                        oGraz  = 12, &
                        oZ2N   = 13, &
                        oD2N   = 14, &
                        odmudl = 15, &
                        od2mu  = 16, &
                        od2gdl = 17, &
                        ow_p   = 18
 
  real(8) :: pp(4,4), PP_pn,Zmort
  logical, parameter :: do_IRON = .false.
!-----------------------------------------------------------------------
!BOC
    
  ! Retrieve current (local) state variable values.
  NO3    = Varin(iNO3)
  PHY    = Varin(iPHY)
  ZOO    = Varin(iZOO)
  DET    = Varin(iDET)
  CHL    = Varin(iCHL)
  PMUPHY = Varin(iPMU)
  VARPHY = Varin(iVAR)
  Varout(:) = 0.
!All rates have been multiplied by dtdays to get the real rate correponding to the actual time step
  tf_p    = temp(Ep,tC)
  DET     = max(DET,eps)
  PHY     = max(PHY,eps)
  PMUPHY  = max(PMUPHY,eps)
  VARPHY  = max(VARPHY,eps)
  PMU     = PMUPHY/PHY
  VAR     = VARPHY/PHY
! Fe related:
  if (do_IRON .eq. .true.) then
     Fe   = Varin(iFER)
     Fe   = max(Fe,Femin)  !Dissolved Fe concentration
     KFe  = ScaleTrait(PMU, K0Fe,alphaFe) !The half saturation constant of Fe at average size
  endif
    
  Qs      = ScaleTrait(PMU, Q0N, alphaQ)/2.
  mu0hat  = dtdays*tf_p*mu0*exp(alphamu*PMU + betamu*PMU*PMU)
  X       = 0.
  
  if (do_IRON .eq. .true.) then
    mu0hat  = mu0hat * Fe/(Fe + KFe)
    X       = alphaFe*KFe/(Fe + KFe) 
  endif
  
  dmu0hatdl = mu0hat*(alphamu + 2D0*betamu*PMU - X)
  
  d2mu0hatdl2= dmu0hatdl*(alphamu + 2D0*betamu*PMU-X) + mu0hat*2D0*betamu
  
  if (do_IRON .eq. .true.) then
     d2mu0hatdl2 = d2mu0hatdl2 - mu0hat*alphaFe*Fe*X/(Fe+KFe)
  endif
  
  ! Iron limits nitrogen uptake
  V0hat  = ScaleTrait(PMU, dtdays*tf_p*V0N, alphaV)
  
  if (do_IRON .eq. .true.) V0hat = V0hat*Fe/(Fe + KFe)
  
  dV0hatdl   = V0hat*(alphaV- X)
  d2V0hatdl2 = dV0hatdl*(alphaV - X)
  
  if (do_Iron .eq. .true.) then
  
    d2V0hatdl2 = d2V0hatdl2 - V0hat*alphaFe*Fe*X/(Fe+KFe) 
  
  endif
  
  ! Initial slope of P-I curve
  aI     = ScaleTrait(PMU, dtdays*aI0, alphaI)
  ! Cost of photosynthesis
  RMchl  = tf_p*RMchl0*dtdays
  ! Threshold irradiance and RMchl is set temperature dependent
  I_zero = zetaChl*RMchl/aI  
  
  !Define VNhat: the saturation function of ambient nutrient concentration
  selectcase(nutrient_uptake)  
  ! case 1: Classic Michaelis Menton 
    case(1)
  ! Potential maximal nutrient-limited uptake
      ! Half-saturation constant of nitrate uptake
      Kn     = ScaleTrait(PMU,K0N, alphaK) 
      VNhat  = V0hat*NO3/(NO3 + Kn)
  
   dVNhatdl  = -VNhat*alphaK*Kn/(NO3+Kn) + NO3/(NO3 + Kn)*dV0hatdl
  
     d2VNhatdl2 = -alphaK*(VNhat * alphaK * NO3/(NO3+Kn)**2*Kn &
   + Kn/(NO3 + Kn)*dVNhatdl)+ NO3/(NO3+ Kn)*d2V0hatdl2         &
   - dV0hatdl*NO3/(NO3 + Kn)**2*alphaK*Kn
  
     ! case 2: optimal uptake based on Pahlow (2005) and Smith et al. (2009)
     case(2)
  
       Lmin  = log(Lref**3/6.*pi) + log(1.+Penfac)/( Penfac*alphaA) 
       A0N   = dtdays*tf_p*A0N0
       A0hat = PenAff(PMU, alphaA, Penfac,Lmin) * ScaleTrait(PMU, A0N, alphaA)
       A0hat = max(A0hat,eps)
  
      dA0hatdl = alphaA*A0hat-A0N*exp(PMU*alphaA)*Penfac*alphaA   &
              * exp(Penfac*alphaA *(PMU-Lmin))
    
      d2A0hatdl2 = alphaA*dA0hatdl   &
              - Penfac*alphaA*exp(alphaA*((1.+Penfac)*PMU-Penfac*Lmin)) &
              * (dA0hatdl + A0hat*alphaA*(1.+Penfac))  
       
              !Define fA
       fA = 1D0/( 1D0 + sqrt(A0hat * NO3/V0hat) ) 
    
       VNhat = (1D0-fA)*V0hat*fA*A0hat*NO3/((1D0-fA)*V0hat + fA*A0hat*NO3) 
       
       !X: temporary variable
       X    = V0hat/A0hat + 2.*sqrt(V0hat*NO3/A0hat) + NO3       
    
       !B: d(V0/A0)dl
       B    = dV0hatdl/A0hat - V0hat/A0hat**2*dA0hatdl
    
       dXdl = B*(1. + sqrt(NO3 * A0hat / V0hat))
    
       dBdl = d2V0hatdl2/A0hat - dV0hatdl*dA0hatdl/A0hat**2    &
        - (V0hat/A0hat**2*d2A0hatdl2      &
        +  dA0hatdl*(dV0hatdl/A0hat**2 - 2.*V0hat*dA0hatdl/A0hat**3))
       
       dVNhatdl = NO3*(dV0hatdl/X-V0hat/X**2*B*(1.+ sqrt(NO3*A0hat/V0hat) ) )
    
       d2VNhatdl2 = NO3*(d2V0hatdl2/X - dV0hatdl*dXdl/X**2    &
        - (V0hat/X**2*B*(-sqrt(NO3)*.5*(A0hat/V0hat)          &
        * sqrt(A0hat/V0hat) * B)                              &
        + B*(1. + sqrt(NO3*A0hat/V0hat) )                     &
        * (dV0hatdl/X**2 - 2.*V0hat*dXdl/X**3)                &
        + V0hat/X**2*(1.+sqrt(NO3*A0hat/V0hat))*dBdl))
    
       case default
        write(*,*) 'Error: Incorrect option for nutrient uptake!'
       ENDSELECT  

! Calculate thetahat (optimal g Chl/mol C for the chloroplast under nutrient replete conditions)
! Only calculate within the euphotic zone, otherwise many numerical problems.

      if( par .gt. I_zero ) then
        
        larg1 = exp(1. + min(aI*par/(mu0hat*zetaChl),600.))
   
        larg  = (1. + RMchl/mu0hat)*larg1   
        
        dmu0hat_aIdl   = (dmu0hatdl - alphaI*mu0hat)/aI
   
        d2mu0hat_aIdl2 = d2mu0hatdl2/aI -alphaI/aI*dmu0hatdl-aI*dmu0hat_aIdl
   
        daI_mu0hatdl = -(aI/mu0hat)**2*dmu0hat_aIdl
   
        d2aI_mu0hatdl2 = -((aI/mu0hat)**2*d2mu0hat_aIdl2   &
       - 2./(mu0hat/aI)**3*(dmu0hat_aIdl)**2)
   
        dlargdl = -RMchl*larg1/(mu0hat**2)*dmu0hatdl       &
       + (1.+RMchl/mu0hat)*larg1 * par/zetaChl*daI_mu0hatdl
        
        d2largdl2 = -RMchl*(larg1*mu0hat**(-2)*d2mu0hatdl2   &
       + larg1*par/zetaChl*daI_mu0hatdl*mu0hat**(-2)*dmu0hatdl  &
       + larg1*dmu0hatdl*(-2.*mu0hat**(-3)*dmu0hatdl))   &
       + par/zetaChl*((1+RMchl/mu0hat)*larg1*d2aI_mu0hatdl2  &
       + (1.+RMchl/mu0hat)*larg1*par/zetaChl*daI_mu0hatdl*daI_mu0hatdl &
       + RMchl*(-mu0hat**(-2)*dmu0hatdl)*larg1*daI_mu0hatdl)
   
       W_Y      = WAPR(larg,0,0)
       ThetaHat = 1.0/zetaChl + (1.0- W_Y)*mu0hat/(aI * par)
       ThetaHat = max(ThetaHat,thetamin) 
   
       dThetaHatdl = 1./par*(-W_Y/larg/(1.+W_Y)*dlargdl*mu0hat/aI  &
     +  (1.-W_Y)*dmu0hat_aIdl)
       
       dWYYdl = dlargdl*(-W_Y**2/larg**2/(1.+W_Y)**3   &
       -  W_Y/larg**2/(1.+W_Y) + W_Y/(larg*(1.+W_Y))**2)
   
       d2ThetaHatdl2 = 1/par*(-(W_Y/larg/(1.+W_Y)*dlargdl*dmu0hat_aIdl  &
       +  W_Y/larg/(1.+W_Y)*d2largdl2*mu0hat/aI   &
       +  dWYYdl*dlargdl*mu0hat/aI)               &
       -  W_Y/larg/(1.+W_Y)*dlargdl * dmu0hat_aIdl &
       +  (1.-W_Y)*d2mu0hat_aIdl2)
   
        SI = 1. - max(exp(-aI*par*ThetaHat/mu0hat),0.)
   
        dSIdl = ( (alphaI- dmu0hatdl/mu0hat)   &
        * ThetaHat + dThetaHatdl) * (1.-SI)*aI*par/mu0hat    !confirmed
   
       d2SIdl2 = par*(- dSIdl*aI/mu0hat*(ThetaHat*alphaI     &
      - ThetaHat/mu0hat*dmu0hatdl + dThetaHatdl) + (1.-SI)   &
      * (ThetaHat*alphaI- ThetaHat/mu0hat*dmu0hatdl + dThetaHatdl) &
      * daI_mu0hatdl + (1.-SI)*aI/mu0hat*(                 &
      - (d2mu0hatdl2/mu0hat - dmu0hatdl**2/mu0hat**2)*ThetaHat  &
      + (alphaI-dmu0hatdl/mu0hat)*dThetaHatdl + d2ThetaHatdl2)  )
   
       ! Light dependent growth rate 
       ! (needs to take into account the cost of dark and light-dependent chl maintenance)
       mu0hatSI = mu0hat*SI  ! Gross specific carbon uptake (photosynthesis)
   
       muIhat   = mu0hatSI-(mu0hatSI+RMchl)*zetaChl*ThetaHat ! Net specific carbon uptake
       muIhat   = max(muIhat,eps)
   
       dmu0hatSIdl = SI*dmu0hatdl + mu0hat*dSIdl
   
       d2mu0hatSIdl2 = d2mu0hatdl2*SI+2.*dmu0hatdl*dSIdl+mu0hat*d2SIdl2 !Correct
   
       dmuIhatdl = (1D0-zetaChl*ThetaHat)*dmu0hatSIdl - dThetaHatdl*zetaChl*(mu0hatSI+RMchl) !Correct
    
       d2muIhatdl2=d2mu0hatSIdl2 - zetaChl*(ThetaHat*d2mu0hatSIdl2+2.*dThetaHatdl*dmu0hatSIdl   &
       +mu0hatSI*d2ThetaHatdl2)-zetaChl*RMchl*d2ThetaHatdl2   !Correct
    
       ZINT   = Qs*(muIhat/VNhat + zetaN)
    
       dZINdl = Qs*(dmuIhatdl/VNhat - muIhat*dVNhatdl/VNhat**2)+alphaQ*ZINT    
       
       d2ZINdl2 = Qs/VNhat*((alphaQ-dVNhatdl/VNhat)*dmuIhatdl & 
       + d2muIhatdl2) - Qs/VNhat**2*(muIhat*d2VNhatdl2        &
       + dVNhatdl*(dmuIhatdl+alphaQ*muIhat-2.*muIhat          &
       / VNhat*dVNhatdl)) + alphaQ*dZINdl
       
       fV = (-1.0 + sqrt(1.0 + 1.0/ZINT))*Qs*muIhat/VNhat
       fV = max(fV,fVmin)
    !
       else
    ! Under the conditions of no light:
          ThetaHat      = thetamin  !  a small positive value 
          dThetaHatdl   = 0.
          d2ThetaHatdl2 = 0.
          ZINT          = Qs*zetaN
          dZINdl        = alphaQ*ZINT
          d2ZINdl2      = alphaQ*dZINdl
          fV            = fVmin
          muIhat        = 0.
          dmuIhatdl     = 0.
          d2muIhatdl2   = 0.
       endif
    
           ! Optimal nutrient quota:
        QN = (1.+ sqrt(1.+1./ZINT))*Qs
    
        dQNdl  = alphaQ*QN-dZINdl*Qs/(2.*ZINT*sqrt(ZINT*(1.+ZINT))) !confirmed  
    
        d2QNdl2 = alphaQ*dQNdl - Qs/(2.*ZINT*sqrt(ZINT*(ZINT+1.)))  &
         *(d2ZINdl2+alphaQ*dZINdl-(2.*ZINT+1.5)/(ZINT*(ZINT+1.))*dZINdl**2)      ! Confirmed
    
        dfVdl = alphaQ*Qs*(1/QN+2.*zetaN)-(zetaN+Qs/QN**2)*dQNdl  !Confirmed
    !
        d2fVdl2 = (alphaQ**2)*Qs*(1/QN + 2*zetaN)                    &
          -  2.*alphaQ*Qs*dQNdl/QN**2 +  2.*(dQNdl**2)*Qs/QN**3      &     
          -     (zetaN + Qs/QN**2) * d2QNdl2  ! Confirmed
    
    
            if (par .gt. I_zero) then
    !Net growth rate (d-1) of phytoplankton at the average size
                muNet = muIhat*(1.-fV-Qs/QN) - zetaN*fV*VNhat
    !Here the derivative of muNet includes respiratory costs of both N Assim and Chl maintenance       
              dmuNetdl = muIhat*(Qs/(QN**2)*dQNdl                         &
          -  alphaQ*Qs/QN-dfVdl) +  (1.-fV-Qs/QN)*dmuIhatdl  & 
          -  zetaN*(fV*dVNhatdl+VNhat*dfVdl)
    
             d2muNetdl2 = (Qs/(QN*QN))*dQNdl*dmuIhatdl    &
         + muIhat*(Qs/QN**2)*(dQNdl*(alphaQ - 2.*dQNdl/QN) + d2QNdl2)  & 
         - (alphaQ*(Qs/QN*(dmuIhatdl + alphaQ*muIhat)  &
         - muIhat*Qs/(QN**2)*dQNdl))                   &
         - (muIhat*d2fVdl2+dfVdl*dmuIhatdl)            &
         + dmuIhatdl*(Qs/(QN**2)*dQNdl-alphaQ*Qs/QN-dfVdl)  & 
         + (1.-fV-Qs/QN)*d2muIhatdl2                   &
         - zetaN*(fV*d2VNhatdl2 + 2.*dfVdl*dVNhatdl + VNhat*d2fVdl2)  !dC
    
            else
             muNet      = 0.
             dmuNetdl   = 0.
             d2muNetdl2 = 0.
            endif
    !  chl:C ratio [ g chl / mol C ] of the whole cell, at the mean cell size 
       Theta = ThetaHat*(1. - fV - Qs/QN)
    ! Calculate the mean chl:C ratio of phytoplankton (averaged over all sizes) 
    ! using the second derivative of the cell quota w.r.t. log size. (Why negative?)
       ThetaAvg = max(Theta + VAR*.5*(d2ThetaHatdl2*(1.0 - fV - Qs/QN)  &     
       - ThetaHat*(d2fVdl2 + Qs*(dQNdl**2*(2./QN)-d2QNdl2)/(QN**2))),thetamin)
                
    ! Calculate the mean N:C ratio of phytoplankton (averaged over all sizes) 
    ! using the second derivative of the cell quota w.r.t. log size. (How to derive?)
       QNAvg = max(QN/(1.+((2/QN)*dQNdl**2 - d2QNdl2)*VAR/(2.*QN)),QNmin)  
    
!=============================================================
   ! Calculate the community sinking rate of phytoplankton
   ! Phytoplankton sinking rate at the average size
   w_p    = ScaleTrait(PMU,abs(dtdays*w_p0),alphaW)  !Positive
   ! Constrain the sinking rate not too large for large cells (Smayda 1970)
   w_p    = min(w_p, Wmax*dtdays*1.)
   dwdl   = alphaW*w_p
   d2wdl2 = alphaW*dwdl
 ! Sinking rate (m) of phytoplankton
   w_pAvg = w_p+0.5*VAR*d2wdl2
 ! sinking rate must be negative!
   w_pAvg = -min(w_pAvg, Wmax*dtdays*1.)
    
!=============================================================

!! ZOOplankton section:
   tf_z = temp(Ez,tC)
  
   ! The grazing dependence on total prey
   gbar = grazing(grazing_formulation,kp,PHY)

   !Zooplankton total ingestion rate
   INGES = tf_z*dtdays*gmax*gbar
   !Zooplankton excretion rate (-> DOM)
   RES = INGES*(1.-GGE-unass)
   !ZOOPLANKTON EGESTION (-> POM)
   EGES = INGES*unass
    
! Grazing rate on the mean size of PHY (specific to N-based Phy biomass, unit: d-1) (Eq. 12)
 gbar = INGES*ZOO/PHY*sqrt(alphaG)
 dgdlbar   = 0.
 d2gdl2bar = (1.-alphaG)/VAR*gbar
!!End of zooplankton section
    
!=============================================================
!! Solve ODE functions:
   Zmort = ZOO*ZOO*dtdays*mz*tf_z  !Mortality term for ZOO

   if (d2gdl2bar .lt. 0.) then
      VAR1 = VAR*(1.-VAR*d2gdl2bar)
   else                                        
      VAR1 = VAR/(1.+VAR*d2gdl2bar)
   endif

   ! Update PMU and VAR:    
   if (d2muNetdl2 .gt. 0.) then
! self%dVARdt=VAR*VAR*(self%d2muNetdl2- self%d2gdl2bar-self%d2wdl2)  !Eq. 19 
      VAR1 =VAR1*(1.+VAR*d2muNetdl2)
   else    
      VAR1 =VAR1/(1.-VAR*d2muNetdl2)
   endif
   VAR1    =VAR1/(1.+VAR*d2wdl2)
   
   !Eq. 18, Update PMU:
   !   self%dPMUdt = self%p%VAR*(self%dmuNetdl - self%dgdlbar - self%dwdl)
   if (dmuNetdl .gt. 0.) then
      PMU1 = PMU  + VAR * dmuNetdl
   else
      PMU1 = PMU/(1.-VAR/PMU*dmuNetdl)
   endif
   
   PMU1 = PMU1/(1.+ VAR/PMU*dwdl)
   !Contrain the PMU and VAR: 
   PMU1 = min(max(PMU1,0.01),PMUmax)
   VAR1 = min(max(VAR1,0.01),VARmax)
          
    ! For production/destruction matrix:
      
   pp(:,:) = 0.
   pp(iNO3,iDET)=dtdays*rDN*DET*tf_z   
   pp(iNO3,iZOO)=ZOO*RES        
   pp(iDET,iZOO)=ZOO*EGES+Zmort 
   pp(iZOO,iPHY)=ZOO*INGES      
   pp(iDET,iPHY)=PHY*PHY*tf_p*dtdays*mp
   
   PP_pn = PHY*(muNet + 0.5*VAR*d2muNetdl2)
   pp(iPHY,iNO3)=max(0.5*( PP_pn + abs(PP_pn)), 0.)
   pp(iNO3,iPHY)=max(0.5*(-PP_pn + abs(PP_pn)), 0.)
   
   DET1 = (DET + pp(iDET,iZOO) + pp(iDET,iPHY))/(1.+ pp(iNO3,iDET)/DET )
   DET1 = min(max(DET1,eps),DETmax)
   
   NO31 = (NO3+pp(iNO3,iDET)+pp(iNO3,iZOO)+pp(iNO3,iPHY))/(1.+pp(iPHY,iNO3)/NO3)
    
   PHY1=(PHY+pp(iPHY,iNO3))  &
       /(1.+(pp(iZOO,iPHY)+pp(iNO3,iPHY)+pp(iDET,iPHY))/PHY)
    
   ZOO1 = (ZOO+pp(iZOO,iPHY))/(1.+ EGES+ZOO*dtdays*mz*tf_z+RES)
   
   CHL1 = PHY1/QNavg*ThetaAvg
   PMUPHY1 =(PMUPHY+PHY*PMU1+PMU*PHY1)/(1.+ 2.*PHY*PMU/PMUPHY)
   VARPHY1 =(VARPHY+PHY*VAR1+VAR*PHY1)/(1.+ 2.*PHY*VAR/VARPHY)
   
   Varout(oNO3)   = NO31
   Varout(oPHY)   = PHY1
   Varout(oZOO)   = ZOO1
   Varout(oDET)   = DET1
   Varout(oPMU)   = PMUPHY1
   Varout(oVAR)   = VARPHY1
   Varout(oCHL)   = CHL1

   if(do_IRON .eq. .true.) Varout(oFER)   = Fe
   Varout(oTheta) = ThetaAvg
   Varout(oQN)    = QNAvg
   Varout(omuNet) = PP(iPHY,iNO3)/dtdays/PHY
   Varout(oGraz ) = PP(iZOO,iPHY)/dtdays/PHY
   Varout(oZ2N  ) = PP(iNO3,iZOO)/dtdays
   Varout(oD2N  ) = PP(iNO3,iDET)/dtdays
   Varout(odmudl) = dmuNetdl/dtdays
   Varout(od2mu ) = d2muNetdl2/dtdays
   Varout(od2gdl) = d2gdl2bar/dtdays
   Varout(ow_p)   = w_pAvg/dtdays
   return
END SUBROUTINE FlexEFT 
!========================================================
real(8) function temp(Ea,tC)
implicit none
!DESCRIPTION:
!The temperature dependence of plankton rates are fomulated according to the Arrhenuis equation. 
! tC: in situ temperature
! Tr: reference temperature
!
!INPUT PARAMETERS:
 real(8), intent (in) :: Ea, tC
 ! boltzman constant constant [ eV /K ]
 real(8), parameter   :: kb = 8.62d-5, Tr = 15.

temp = exp(-(Ea/kb)*(1D0/(273.15 + tC)-1D0/(273.15 + Tr)))

return 
end function temp
!========================================================
  FUNCTION WAPR (X, NB, L) RESULT (WAP)
!
!     WAPR - output
!     X - argument of W(X)
!     NB is the branch of the W function needed:
!        NB = 0 - upper branch
!        NB <> 0 - lower branch
!
!     NERROR is the output error flag:
!        NERROR = 0 -> routine completed successfully
!        NERROR = 1 -> X is out of range
!
!     Range: -exp(-1) <= X for the upper branch of the W function
!            -exp(-1) < X < 0 for the lower branch of the W function
!
!     L - determines how WAPR is to treat the argument X
!        L = 1 -> X is the offset from -exp(-1), so compute
!                 W(X-exp(-1))
!        L <> 1 -> X is the desired X, so compute W(X)
!
!     M - print messages from WAPR?
!         M = 1 -> Yes
!         M <> 1 -> No
!
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     NN is the output device number
!
!     NBITS is the number of bits (less 1) in the mantissa of the
!        floating point number number representation of your machine.
!        It is used to determine the level of accuracy to which the W
!        function should be calculated.
!
!        Most machines use a 24-bit matissa for single precision and
!        53-56 bits for double precision. The IEEE standard is 53
!        bits. The Fujitsu VP2200 uses 56 bits. Long word length
!        machines vary, e.g., the Cray X/MP has a 48-bit mantissa for
!        single precision.
!
    IMPLICIT NONE
    real, INTENT(in)   :: X 
    INTEGER, PARAMETER :: NN=6, NBITS=23, NITER=1
    real, PARAMETER ::EM=-0.367879441171442,&           ! -EXP(-1)
                      EM9=-1.234098040866796E-4,&       ! -EXP(-9)
                      C13=1.E0/3.E0,&
                      C23=2.E0*C13,&
                      EM2=2.E0/EM,&
                      E12=-EM2,&
                      TB=.5E0**NBITS,&
                      TB2=.5E0**(NBITS/2),&       ! SQRT(TB)
                      X0=0.0350769390096679055,&  ! TB**(1/6)*0.5E0
                      X1=0.302120119432784731,&   !(1 - 17*TB**(2/7))*EM
                      AN22=3.6E2/83.E0,&
                      AN11=135./83.E0,&
                      AN3=8.E0/3.E0,&
                      AN4=135.E0/83.E0,&
                      AN5=166.E0/39.E0,&
                      AN6=3167.E0/3549.E0,&
                      S2=1.41421356237310,& ! SQRT(2.E0)
                      S21=2.E0*S2-3.E0,&
                      S22=4.E0-3.E0*S2,&
                      S23=S2-2.E0
    real ::  WAP, AN2, DELX,  RETA, ZL, T, TS, ETA, TEMP, TEMP2, ZN
    real ::  XX
    INTEGER, INTENT(in) :: NB, L

!        Various mathematical constants
    
!
!     The following COMMON statement is needed only when testing this
!     function using BISECT, otherwise it can be removed.
!
!    COMMON/WAPCOM/NBITS
!    DATA INIT,NITER/0,1/
!     DATA NBITS/23/
!
!     IF(INIT.EQ.0) THEN
!        INIT=1
!
!        Code to calculate NBITS for the host machine. NBITS is the
!        mantissa length less one. This value is chosen to allow for
!        rounding error in the final bit. This need only be run once on
!        any particular machine. It can then be included in the above
!        DATA statement.
!
!        DO I=1,2000
!           B=2.E0**(-I)
!           V=1.E0+B
!           IF(V.EQ.1.E0)THEN
!              NBITS=I-1
!              J=-ALOG10(B)
!              IF(M.EQ.1) WRITE(NN,40)NBITS,J
!              EXIT
!           ENDIF
!        END DO
!
!        Remove to here after NBITS has been calculated once
!
!        The case of large NBITS
!
!        IF(NBITS.GE.56) NITER=2
!
!        Various mathematical constants
!
!        EM=-EXP(-1.E0)
!        EM9=-EXP(-9.E0)
!        C13=1.E0/3.E0
!        C23=2.E0*C13
!        EM2=2.E0/EM
!        E12=-EM2
!        TB=.5E0**NBITS
!        TB2=SQRT(TB)
!        X0=TB**(1.E0/6.E0)*.5E0
!        X1=(1.E0-17.E0*TB**(2.E0/7.E0))*EM
!        AN22=3.6E2/83.E0
!        AN11=135./83.E0
!        AN3=8.E0/3.E0
!        AN4=135.E0/83.E0
!        AN5=166.E0/39.E0
!        AN6=3167.E0/3549.E0
!        S2=SQRT(2.E0)
!        S21=2.E0*S2-3.E0
!        S22=4.E0-3.E0*S2
!        S23=S2-2.E0
!     ENDIF
    IF(L.EQ.1) THEN
       DELX=X
       IF(DELX.LT.0.E0) THEN
          WAP = 1./0.
          RETURN
       END IF
       XX=X+EM
!        IF(E12*DELX.LT.TB**2.AND.M.EQ.1) WRITE(NN,60)DELX
    ELSE
       IF(X.LT.EM) THEN
          WAP = 1./0.
          RETURN
       END IF
       IF(X.EQ.EM) THEN
          WAP=-1.E0
          RETURN
       ENDIF
       XX=X
       DELX=XX-EM
!        IF(DELX.LT.TB2.AND.M.EQ.1) WRITE(NN,70)XX
    ENDIF
    IF(NB.EQ.0) THEN
!
!        Calculations for Wp
!
       IF(ABS(XX).LE.X0) THEN
          WAP=XX/(1.E0+XX/(1.E0+XX/(2.E0+XX/(.6E0+.34E0*XX))))
          RETURN
       ELSE IF(XX.LE.X1) THEN
          RETA=SQRT(E12*DELX)
          WAP=RETA/(1.E0+RETA/(3.E0+RETA/(RETA/(AN4+RETA/(RETA*&
               AN6+AN5))+AN3)))-1.E0
          RETURN
       ELSE IF(XX.LE.2.E1) THEN
          RETA=S2*SQRT(1.E0-XX/EM)
          AN2=4.612634277343749E0*SQRT(SQRT(RETA+&
               1.09556884765625E0))
          WAP=RETA/(1.E0+RETA/(3.E0+(S21*AN2+S22)*RETA/&
               (S23*(AN2+RETA))))-1.E0
       ELSE
          ZL =ALOG(XX)
          WAP=ALOG(XX/ALOG(XX/ZL**EXP(-1.124491989777808E0/&
               (.4225028202459761E0+ZL))))
       ENDIF
    ELSE
!
!        Calculations for Wm
!
       IF(XX.GE.0.E0) THEN
          WAP = -1./0.
          RETURN
       END IF
       IF(XX.LE.X1) THEN
          RETA=SQRT(E12*DELX)
          WAP=RETA/(RETA/(3.E0+RETA/(RETA/(AN4+RETA/(RETA*&
               AN6-AN5))-AN3))-1.E0)-1.E0
          RETURN
       ELSE IF(XX.LE.EM9) THEN
          ZL=ALOG(-XX)
          T=-1.E0-ZL
          TS=SQRT(T)
          WAP=ZL-(2.E0*TS)/(S2+(C13-T/(2.7E2+&
               TS*127.0471381349219E0))*TS)
       ELSE
          ZL=ALOG(-XX)
          ETA=2.E0-EM2*XX
          WAP=ALOG(XX/ALOG(-XX/((1.E0-.5043921323068457E0*&
               (ZL+1.E0))*(SQRT(ETA)+ETA/3.E0)+1.E0)))
       ENDIF
    ENDIF
!     DO I=1,NITER
       ZN=ALOG(XX/WAP)-WAP
       TEMP=1.E0+WAP
       TEMP2=TEMP+C23*ZN
       TEMP2=2.E0*TEMP*TEMP2
       WAP=WAP*(1.E0+(ZN/TEMP)*(TEMP2-ZN)/(TEMP2-2.E0*ZN))
!     END DO
    RETURN
  END FUNCTION WAPR

!====================================================
  real(8) function ScaleTrait( logsize, star, alpha ) 
     implicit none
     real(8), intent(IN) :: logsize, star, alpha
  
     ! Calculate the size-scaled value of a trait
     ! for the given log (natural, base e) of cell volume as pi/6*ESD**3 (micrometers). 
  
     ScaleTrait = star * exp( alpha * logsize )
  
     return
  end function ScaleTrait

!====================================================

real(8) function PenAff( logsize, alpha, Pfac, lmin ) 
  implicit none
  real(8), intent(IN) :: logsize, alpha, Pfac, lmin 

!A 'penalty' function to reduce the value of affinity for nutrient at very small cell sizes
!in order to avoid modeling unrealistically small cell sizes.  This is needed because affnity
!increases with decreasing cell size, which means that under low-nutrient conditions, without
!such a penalty, unrealistically small cell sizes could be predicted.
!This penalty function becomes zero at logsize = lmin.   
   
  PenAff = 1.0 - exp(Pfac*alpha*(logsize - lmin))
end function PenAff
!====================================================

real(8) function grazing(Hollingtype, Ksat, Prey)
  implicit none
  real(8),intent(in) :: Ksat, Prey
  integer, intent(in) :: Hollingtype
  ! kp relates to the capture coefficient
  SELECT CASE(Hollingtype)
    ! Holling Type I
    case (1)
      grazing = min(Prey/2.0/Ksat,1.0)
    ! Holling Type II
    case (2)
      grazing = Prey/(Ksat + Prey)  
    ! Holling Type III
    case (3) 
      grazing = Prey*Prey/(Ksat*Ksat + Prey*Prey)
   ! Ivlev
    case (4)
   !To be consistent with other functions  
      grazing = 1.-exp(-log(2d0)*Prey/Ksat) 

  END SELECT
  return
end function
!===========
subroutine Calculate_PAR(I_0, nlev, Hz, Chl, PAR)
  !top level is nlev, bottom layer is 1, following ROMS convention
  implicit none
  real(8), intent(in) :: I_0
  integer, intent(in) :: nlev    ! Total number of vertical layers
  real(8), intent(in) :: Hz(nlev), Chl(nlev)
  real(8), intent(out):: PAR(nlev)
  integer             :: i
  real(8)             :: par0, attn    ! Scratch variable
  real(8), parameter  :: kw = 0.04, kc = 0.025

  par0 = I_0   !Light at the grid surface
  do i = nlev,1,-1
     attn   = exp(-0.5*Hz(i)*(kw + kc*Chl(i)))
     PAR(i) = par0*attn
     par0   = PAR(i)*attn
  enddo

end subroutine
!===========
! !ROUTINE: Advection schemes --- grid centers\label{sec:advectionMean}
subroutine adv_center(N,dt,h,ho,ww,Bcup,Bcdw,Yup,Ydw,method,mode,Y)
   IMPLICIT NONE
! !INPUT PARAMETERS:
!  number of vertical layers
   integer,  intent(in)                :: N
!  time step (s)
   real(8), intent(in)                 :: dt
!  layer thickness (m)
   real(8), intent(in)                 :: h(N)

!  old layer thickness (m)
   real(8), intent(in)                 :: ho(N)

!  vertical advection speed
   real(8), intent(in)                 :: ww(0:N)

!  type of upper BC
   integer,  intent(in)                :: Bcup

!  type of lower BC
   integer,  intent(in)                :: Bcdw

!  value of upper BC
   real(8), intent(in)                 :: Yup

!  value of lower BC
   real(8), intent(in)                 :: Ydw

!  type of advection scheme
   integer,  intent(in)                :: method

!  advection mode (0: non-conservative, 1: conservative)
   integer,  intent(in)                :: mode
!
! !INPUT/OUTPUT PARAMETERS:
   real(8),  intent(inout)             :: Y(N)
!
! !DEFINED PARAMETERS:
   real(8),      parameter             :: one6th=1.0d0/6.0d0
   integer,      parameter             :: itmax=100

!  type of advection scheme
   integer,parameter                    :: UPSTREAM       = 1
   integer,parameter                    :: P1             = 2
   integer,parameter                    :: P2             = 3
   integer,parameter                    :: Superbee       = 4
   integer,parameter                    :: MUSCL          = 5
   integer,parameter                    :: P2_PDM         = 6

!  boundary condition type
!  for advection schemes
   integer,parameter                    :: flux           = 1
   integer,parameter                    :: value          = 2
   integer,parameter                    :: oneSided       = 3
   integer,parameter                    :: zeroDivergence = 4

! !LOCAL VARIABLES:
   integer                             :: i,k,it
   real(8)                             :: x,r,Phi,limit
   real(8)                             :: Yu,Yc,Yd
   real(8)                             :: c,cmax
   real(8)                             :: cu(0:N)
!-----------------------------------------------------------------------
!  initialize interface fluxes with zero
   cu   = 0.0

!  initialize maximum Courant number
   cmax = 0.0

!  compute maximum Courant number
   do k = 1,N-1
      c = abs(ww(k))*dt/(0.5*(h(k)+h(k+1)))
      if (c .gt. cmax) cmax=c
   enddo

   it   = min(itmax,int(cmax)+1)

!#ifdef DEBUG
   if (it .gt. 1) then
      write(0,*) 'In adv_center():'
      write(0,*) 'Maximum Courant number is ',cmax
      write(0,*) it,' iterations used for vertical advection'
   endif
!#endif

!  splitting loop
   do i=1,it

!     vertical loop
      do k=1,N-1

!        compute the slope ration
         if (ww(k) .gt. 0.0) then

!           compute Courant number
            c=ww(k)/float(it)*dt/(0.5*(h(k)+h(k+1)))

            if (k .gt. 1) then
               Yu=Y(k-1)                              ! upstream value
            else
               Yu=Y(k)
            end if
            Yc=Y(k  )                                 ! central value
            Yd=Y(k+1)                                 ! downstream value

!           compute slope ratio
            if (abs(Yd-Yc) .gt. 1e-10) then
               r=(Yc-Yu)/(Yd-Yc)
            else
               r=(Yc-Yu)*1.e10
            end if

!        negative speed
         else

!           compute Courant number
            c=-ww(k)/float(it)*dt/(0.5*(h(k)+h(k+1)))

            if (k .lt. N-1) then
               Yu=Y(k+2)                              ! upstream value
            else
               Yu=Y(k+1)
            end if
            Yc=Y(k+1)                                 ! central value
            Yd=Y(k  )                                 ! downstream value

!           compute slope ratio
            if (abs(Yc-Yd) .gt. 1e-10) then
               r=(Yu-Yc)/(Yc-Yd)
            else
               r=(Yu-Yc)*1.e10
            end if

         end if

!        compute the flux-factor phi
         x    =  one6th*(1.-2.0*c)
         Phi  =  (0.5+x)+(0.5-x)*r

!        limit the flux according to different suggestions (Pietrzak 1998)
         select case (method)
            case (UPSTREAM)
               limit= 0.0
            case (P1)
               write(0,*) "Fatal error: P1 advection method not yet implemented, choose other method"
               stop  "adv_center.F90"
            case ((P2),(P2_PDM))
               if (method .eq. P2) then
                  limit=Phi
               else
                  limit=max(0.0,min(Phi,2./(1.-c),2.*r/(c+1.e-10)))
               end if
            case (Superbee)
               limit=max(0.0, min(1d0, 2.0*r), min(r,2.*1d0) )
            case (MUSCL)
               limit=max(0.0,min(2.*1d0,2.0*r,0.5*(1.0+r)))
            case default
               write(0,*) method
               write(0,*) 'Fatal error: unkown advection method in adv_center()'
               stop
          end select

!        compute the limited flux
         cu(k)=ww(k)*(Yc+0.5*limit*(1-c)*(Yd-Yc))

      end do

!     do the upper boundary conditions
      select case (Bcup)
      case (flux)
         cu(N) = - Yup              ! flux into the domain is positive
      case (value)
         cu(N) =  ww(N)*Yup
      case (oneSided)
         if (ww(N).ge.0.0) then
            cu(N) =  ww(N)*Y(N)
         else
            cu(N) = 0.0
         end if
      case (zeroDivergence)
         cu(N) = cu(N-1)
      case default
         write(0,*) 'Fatal error: unkown upper boundary condition type in adv_center()'
         stop
      end select

!     do the lower boundary conditions
      select case (Bcdw)
      case (flux)
         cu(0) =   Ydw               ! flux into the domain is positive
      case (value)
         cu(0) =  ww(0)*Ydw
      case (oneSided)
         if(ww(0).le.0.0) then
            cu(0) =  ww(0)*Y(1)
         else
            cu(0) = 0.0
         end if
      case (zeroDivergence)
         cu(0) = cu(1)
      case default
         write(0,*) 'Fatal error: unkown lower boundary condition type in adv_center()'
         stop
      end select

!     do the vertical advection step which will be used for prescribed
!     vertical flow velocity and for settling of suspended matter.

      if (mode.eq.0) then ! non-conservative, water vertical velocity
         do k=1,N
            Y(k)=Y(k)-1./float(it)*dt*((cu(k)-cu(k-1))/        &
                 h(k)-Y(k)*(ww(k)-ww(k-1))/h(k))
         enddo
      else                ! conservative, for PHY and DET sinking
         do k=1,N
            Y(k)=Y(k)-1./float(it)*dt*((cu(k)-cu(k-1))/h(k))
         enddo
      end if

   end do ! end of the iteration loop
  return
end subroutine adv_center
!EOC
!-----------------------------------------------------------------------
subroutine diff_center(N,dt,cnpar,posconc,h,Bcup,Bcdw, &
                       Yup,Ydw,nuY,Lsour,Qsour,Taur,Yobs,Yin,Yout)

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,parameter                    :: Dirichlet      = 0
   integer,parameter                    :: Neumann        = 1

!  number of vertical layers
   integer,  intent(in)               :: N

!  time step (s)
   real(8), intent(in)                :: dt

!  "implicitness" parameter
   real(8), intent(in)                :: cnpar

!  1: non-negative concentration, 0: else
   integer, intent(in)                :: posconc

!  layer thickness (m)
   real(8), intent(in)                :: h(N)

!  type of upper BC
   integer,  intent(in)               :: Bcup

!  type of lower BC
   integer,  intent(in)               :: Bcdw

!  value of upper BC
   real(8), intent(in)                :: Yup

!  value of lower BC
   real(8), intent(in)                :: Ydw

!  diffusivity of Y
   real(8), intent(in)                :: nuY(0:N)

!  linear source term
!  (treated implicitly)
   real(8), intent(in)                :: Lsour(N)

!  constant source term
!  (treated explicitly)
   real(8), intent(in)                :: Qsour(N)

!  relaxation time (s)
   real(8), intent(in)                :: Taur(N)

!  observed value of Y
   real(8), intent(in)                :: Yobs(N)
!
! !INPUT/OUTPUT PARAMETERS:
   real(8), intent(in)                :: Yin(N)
   real(8), intent(out)               :: Yout(N)

   real(8), dimension(N)              :: au,bu,cu,du
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
   real(8)                   :: a,c,l
!
!-----------------------------------------------------------------------
!BOC
!
!  Initialize au, bu, cu, du
   au(:) = 0.; bu(:) = 0.; cu(:) = 0.; du(:) = 0.
!  set up matrix
   do i=2,N-1
      c     = 2.0*dt*nuY(i)  /(h(i)+h(i+1))/h(i)
      a     = 2.0*dt*nuY(i-1)/(h(i)+h(i-1))/h(i)
      l     =     dt*Lsour(i)

      cu(i) =-cnpar*c
      au(i) =-cnpar*a
      bu(i) = 1d0 + cnpar*(a + c) - l
      du(i) = (1d0 - (1d0-cnpar)*(a + c))*Yin(i)                  &
            + (1d0 - cnpar)*( a*Yin(i-1) + c*Yin(i+1) ) + dt*Qsour(i)
   end do

!   set up upper boundary condition
   select case(Bcup)
   case(Neumann)
      a     = 2.0*dt*nuY(N-1)/(h(N)+h(N-1))/h(N)
      l     = dt*Lsour(N)

      au(N) =-cnpar*a
      if (posconc .eq. 1 .and. Yup.lt.0d0) then ! Patankar (1980) trick
         bu(N) =  1d0 - au(N) - l  - dt*Yup/Yin(N)/h(N)
         du(N) = Yin(N) + dt*Qsour(N)   &
               + (1d0 - cnpar)*a*(Yin(N-1)-Yin(N))
      else
         bu(N) =  1d0 - au(N) - l
         du(N) = Yin(N) + dt*(Qsour(N)+Yup/h(N))   &
               + (1d0 - cnpar)*a*(Yin(N-1)-Yin(N))
      end if
   case(Dirichlet)
      au(N) = 0d0
      bu(N) = 1d0
      du(N) = Yup
   case default
      write(6,*) 'Fatal error: invalid boundary condition type for upper boundary'
      stop  'diff_center.F90'
   end select

!   set up lower boundary condition
   select case(Bcdw)
   case(Neumann)
      c     = 2.0*dt*nuY(1)/(h(1)+h(2))/h(1)
      l     = dt*Lsour(1)

      cu(1) =-cnpar*c
      if (posconc.eq.1 .and. Ydw.lt.0d0) then ! Patankar (1980) trick
         bu(1) = 1d0 - cu(1) - l - dt*Ydw/Yin(1)/h(1)
         du(1) = Yin(1) + dt*(Qsour(1))   &
               + (1d0 - cnpar)*c*(Yin(2)-Yin(1))
      else
         bu(1) = 1d0 - cu(1) - l
         du(1) = Yin(1) + dt*(Qsour(1)+Ydw/h(1))   &
               + (1d0 - cnpar)*c*(Yin(2)-Yin(1))
      end if
   case(Dirichlet)
      cu(1) = 0d0
      bu(1) = 1d0
      du(1) = Ydw
   case default
      write(6,*) 'Fatal error: invalid boundary condition type for lower boundary'
      stop  'diff_center.F90'
   end select


!  relaxation to observed value
   if (minval(Taur).lt.1.E10) then
      do i=1,N
         bu(i)=bu(i)+dt/Taur(i)
         du(i)=du(i)+dt/Taur(i)*Yobs(i)
      end do
   end if


!  solve linear system
   call tridiagonal(N,au,bu,cu,du,1,N,Yout)

   return
end subroutine diff_center
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine tridiagonal(N,au,bu,cu,du,fi,lt,value)
!
! !DESCRIPTION:
! A linear equation with tridiagonal matrix structure is solved here. The main
! diagonal is stored on {\tt bu}, the upper diagonal on {\tt au}, and the
! lower diagonal on {\tt cu}, the right hand side is stored on {\tt du}.
! The method used here is the simplified Gauss elimination, also called
! \emph{Thomas algorithm}.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: N,fi,lt
   real(8), dimension(N), intent(in)   :: au, bu, cu, du
!
! !OUTPUT PARAMETERS:
   real(8), intent(out)                :: value(N)
!
!EOP
!
! !LOCAL VARIABLES:
   integer                             :: i
   real(8), dimension(N)               :: ru, qu
!
!-----------------------------------------------------------------------
!BOC
   ru(lt)=au(lt)/bu(lt)
   qu(lt)=du(lt)/bu(lt)

   do i=lt-1,fi+1,-1
      ru(i)=au(i)/(bu(i)-cu(i)*ru(i+1))
      qu(i)=(du(i)-cu(i)*qu(i+1))/(bu(i)-cu(i)*ru(i+1))
   end do

   qu(fi)=(du(fi)-cu(fi)*qu(fi+1))/(bu(fi)-cu(fi)*ru(fi+1))

   value(fi)=qu(fi)

   do i=fi+1,lt
      value(i)=qu(i)-ru(i)*value(i-1)
   end do


   return
end subroutine tridiagonal
!=======================================================
subroutine  Readcsv(filename,nrow,ncol,dat)
  implicit none
  character(LEN=20),intent(in) :: filename
  character(LEN=10),dimension(ncol) :: header
  integer,          intent(in) :: nrow, ncol
  real(8), dimension(nrow,ncol), intent(out):: dat
  integer    :: err, ix
  integer,parameter :: stdout=6

   dat(:,:) = 0.0
   OPEN(unit=10,file=filename,status='OLD',iostat=err)

   IF (err /= 0) THEN
     write(stdout,*) 'open ', TRIM(filename),' fails'
     stop
     close(10)

   ELSE
     ! write(stdout,*) 'Start reading ',TRIM(filename),' !'
     read(10,*,iostat=err) header  !read the first row
     if (err .ne. 0) then 
       write(stdout,*) 'End of file reached at Row 1.',     &
                      'The error is: ', err
       stop
     else
       continue
     endif

     do ix = 1,nrow  
       read(10,*,iostat=err) dat(ix,:)
       if (err .ne. 0) then 
         write(stdout,*) 'End of file reached at Row: ',ix
         exit
       else
         continue 
       endif
     enddo  

     ENDIF
     CLOSE(unit=10)
! End reading data
 end subroutine Readcsv
!=======================================================
