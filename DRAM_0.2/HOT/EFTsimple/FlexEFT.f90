MODULE FlexEFT1D
  use BIO_rhs
  implicit none

  private

  public FlexEFT1D_main

  character(LEN=2), public, parameter :: Stn = 'S1'

  ! Model options:
  integer,          public, parameter :: Model_ID = EFTdiscrete

contains

subroutine FlexEFT1D_main
  implicit none

!  namelist /model_setup/    Stn,hmax,thetaS,dtsec,Model_ID,maxD
!  namelist /input_data/     Input_dir,N_NO3,N_Temp,N_Aks
!  namelist /output_data/    Output_dir          

  ! Input files
!  character(LEN=22), parameter :: Input_dir = '~/Working/FlexEFT1D/S1' 

  character(LEN=20)  :: Tempfile,PARfile,Aksfile,NO3file

  real(8), parameter :: hmax = 2.5d2, &   ! Total water depth
                       thetaS= 2d0  , &   ! surface stretching parameter
                       dtsec = 3D2  , &   ! time step in seconds
                     d_per_s = 864d2, &   ! how many seconds in one day
                       Yup   = 0d0,   &
                       Ydw   = 0d0 

                        !how many seconds of one year
  integer, parameter :: y_per_s = INT(d_per_s*360), &
                        !Number of vertical points of NO3,Temp, and Aks profile
                        N_NO3   = 14,               &  
                        N_Temp  = 24,               &
                        N_Aks   = 40,               &
                        nsave   = INT(d_per_s)/INT(dtsec), &! Timesteps to save
                        maxD    = 360,              &  ! Number of days to run
                        maxN    = maxD*INT(d_per_s)/INT(dtsec)  ! Number of time steps
                                 
                  
  real(8),parameter  :: cnpar = 6d-1

  real(8) :: Z_r(1:nlev), Z_w(0:nlev), Hz(nlev),  &  ! Grid variables
             obs_NO3(N_NO3  , 2),                 &
             obs_PAR(1      ,13),                 &
             obs_Temp(N_Temp,13),                 &
             obs_Aks(N_Aks  ,13),                 &
             obs_time(12), &                ! observation time indices 
             time, I_0(1), Aks(0:nlev),           &
             NO3(nlev,1),                         &
             VTemp(nlev,12), VAks(0:nlev,12)

  integer :: i, it, k, j, time_day

  character(LEN=20)    :: outfile

  character(LEN=1)     :: YN
  ! New calculated state variables
  real(8), allocatable :: NewVars(:,:), ww(:,:)             
  real(8)              :: t1, t2, cff

  logical, parameter   :: savefile = .FALSE.

!  call CPU_Time(t1)
  
  ! Data input files:
  Tempfile = trim(Stn)//'_Temp.dat'
  PARfile  = trim(Stn)//'_par.dat'
  Aksfile  = trim(Stn)//'_Aks.dat'
  NO3file  = trim(Stn)//'_NO3.dat' 

  !the fraction of a time step in one day
  dtdays   = dtsec/d_per_s

  ! Select biological model
  call choose_model(Model_ID)

  allocate(NewVars(NVAR,nlev))
 ! Sinking rate
  allocate(ww(0:nlev,NVsinkterms))

  ! Setup grid
  call setup_grid(thetaS, nlev, hmax, Z_r, Z_w, Hz)

  ! Initialize state variables
  call Readcsv(NO3file, N_NO3,  2, obs_NO3) 

  ! Interpolate initial NO3:
  call gridinterpol(N_NO3,1,obs_NO3(:,1),obs_NO3(:,2),   &
                    nlev, Z_r, NO3(:,1)) 

  ! Initialize initial NO3:
  Vars(iNO3,:) = NO3(:,1)

  ! Initialize other variables:
  do k = 1,nlev
     do i = 1,NPHY
        Vars(iPHY(i), k) = 1D-3
     enddo
     Vars(iZOO,k) = 1D-2
     Vars(iDET,k) = 1D-2
  enddo

  ! Initialize output data:
  Varout(:,:) = 1D-3
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
  if (savefile .eq. .true.) then
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
  endif


1 format(A10,1x,A7, <nlev  >(2x,F12.5))
2 format(I10,1x,I7, <nlev  >(2x,F12.5))
3 format(A10,1x,A7, <nlev+1>(2x,F12.5))
4 format(I10,1x,I7, <nlev+1>(2x,F12.5))

  ! For each time step, read in external environmental data
  ! 'START TIME STEPPING'

  DO it = 1, maxN+1

    !  write(6,*) 'Continue?' 
    !  read(6,*)  YN

    !  if (YN .eq. 'N') stop
  ! Calculate current timing (zero is starting time):
     time = float(it-1)*dtsec

  ! Interpolate Temp data throughout the water column:
     call time_interp(int(time), 12, nlev, obs_time, VTemp, Temp)

  ! Interpolate temporal PAR data: 
     call time_interp(int(time), 12, 1, obs_time,obs_PAR(1,2:13), I_0)

  ! Convert the unit of par to W m-2:
     cff = I_0(1)/4d-1

  ! Calculate light field:
     call Calculate_PAR(cff, nlev, Hz, Varout(oCHLt,:), PAR)
     
  ! Interpolate Aks data throughout the water column:
     call time_interp(int(time), 12, nlev+1, obs_time, VAks, Aks)


  ! Save data to output files:
  if (mod(it, nsave) .eq. 1) then
  ! Calculate model time in days:
     time_day=int(time/d_per_s)

!     write(6,*) '----------------------------------------'
!     write(6,*) 'Day ', time_day
!     write(6,5) 'Surface NO3 =',Vars(iNO3,nlev)
!     write(6,5) 'Surface ZOO =',Vars(iZOO,nlev)
!     write(6,5) 'Surface DET =',Vars(iDET,nlev)
!     write(6,5) 'Surface CHLA=',Varout(oCHLt,nlev)
!     write(6,5) 'Surface muNet(3) =',Varout(omuNet(3),nlev)/dtdays
!5    format(A15, 1x, F10.5)

  ! Save Temp data into the output file:
     write(9+oTemp, 2) it, time_day, Temp

  ! Save PAR data into the output file:
     write(9+oPAR, 2) it, time_day, PAR

  ! Save Aks data into the output file:
     write(9+oAks, 4) it, time_day, Aks


  ! Save state variables and diagnostics:
     do i=4,Nout+3
        write(9+i, 2) it, time_day, Varout(i-3,:)
     enddo


  endif  !==> End of saving results
  
   ! Biological rhs:
  selectcase(Model_ID)
 
    case(EFTdiscrete)
      call FLEXEFT_DISC
 
    case(EFTcont)
      call FLEXEFT_CONT

    case(Geiderdisc) 
      write(6,*) 'Not coded yet. Quit.'
      stop
    case(Geidercont) 
      write(6,*) 'Not coded yet. Quit.'
      stop

    case(EFTsimple)
      write(6,*) 'Not coded yet. Quit.'
      stop

    case(Geidersimple)
      write(6,*) 'Not coded yet. Quit.'
      stop

    case default
      write(6,*) 'Error in choosing biological models! Quit.'
      stop
  endselect

  ! Pass the new state variables to Vars
  do j = 1,NVAR
     Vars(j,:)=Varout(j,:)
  enddo

  ! Diffusion:

  do j = 1,NVAR

     call diff_center(nlev,dtsec,cnpar,1,Hz,Aks,Vars(j,:),NewVars(j,:))

     ! Save diffusion fluxes (normalized to per day)
     Varout(oD_NO3+j-1,:) = (NewVars(j,:) - Vars(j,:))/dtdays
 
     ! Update the state variables:
     Vars(j,:)=NewVars(j,:)
    
  enddo


  ! Sinking:
  ! Initialize sinking rate:
  ww(0,   :) = 0.0
  ww(nlev,:) = 0.0


  do k = 1,nlev-1
     do i=1,NPHY
        ww(k,i) = -Varout(ow_p(i),k)         !Phytoplankton sinking rate  
     enddo
     ww(k,NPHY+1) = -Params(iwDET)             !Detritus      sinking rate
  enddo

!subroutine adv_center(N,dt,h,ho,ww,Bcup,Bcdw,Yup,Ydw,method,mode,Y)
  do j = 1,NVsinkterms
     call adv_center(nlev,dtdays,Hz,Hz,ww(:,j),  &
                     1,1,Yup,Ydw,6,1,Vars(Windex(j),:))
  enddo


  ! Check whether the values are valid:
  do j = 1,NVAR
     do k=1,nlev
        if( (Vars(j,k) .ne. Vars(j,k)) .OR. (Vars(j,k) .le. 0.) ) then
            write(6,*) 'At time step ',it
            write(6,*) 'The variable ',trim(Labelout(j+3)), ' is invalid at depth ',Z_r(k)
            stop
        endif
     enddo 
  enddo

  ENDDO    ! ==> End of time stepping

  ! Close files:
  if (savefile .eq. .true.) then
     do i = 1, Nout
        close (unit=9+i)
     enddo
  endif

!  ! Release memory:
  deallocate(Vars)
  deallocate(NewVars)
  deallocate(ww)
  deallocate(Varout)
  deallocate(params)
  deallocate(iPHY)
  deallocate(Labelout)
  deallocate(Windex)
  deallocate(oPHY)
  deallocate(oTheta)
  deallocate(oQN)
  deallocate(omuNET)
  deallocate(oGraz)
  deallocate(ow_p)
  deallocate(oSI)
  deallocate(oLno3)
  deallocate(oTheHat)
  deallocate(oD_PHY)

!  call CPU_Time(t2)
!
!  write(6,7) 'Model run finishes. CPU time:  ',t2-t1,' seconds'
7 format(A30,F8.1,A8)

END subroutine
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

end subroutine Calculate_PAR
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
         x    =  one6th*(1d0-2d0*c)
         Phi  =  (5d-1+x)+(5d-1-x)*r

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
                  limit=max(0.0,min(Phi,2d0/(1d0-c),2d0*r/(c+1.e-10)))
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
!subroutine diff_center(N,dt,cnpar,posconc,h,Bcup,Bcdw, &
!                       Yup,Ydw,nuY,Lsour,Qsour,Taur,Yobs,Yin,Yout)
subroutine diff_center(N,dt,cnpar,posconc,h,nuY,Yin,Yout)

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!   integer,parameter                    :: Dirichlet      = 0
!   integer,parameter                    :: Neumann        = 1

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
!   integer,  intent(in)               :: Bcup

!  type of lower BC
!   integer,  intent(in)               :: Bcdw

!  value of upper BC
!   real(8), intent(in)                :: Yup

!  value of lower BC
!   real(8), intent(in)                :: Ydw

!  diffusivity of Y
   real(8), intent(in)                :: nuY(0:N)

!  linear source term
!  (treated implicitly)
!   real(8), intent(in)                :: Lsour(N)

!  constant source term
!  (treated explicitly)
!   real(8), intent(in)                :: Qsour(N)

!  relaxation time (s)
!   real(8), intent(in)                :: Taur(N)

!  observed value of Y
!   real(8), intent(in)                :: Yobs(N)
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
!   Debug:
!   write(6,*) 'Yup    = ',Yup
!   write(6,*) 'Ydw    = ',Ydw

!  Initialize au, bu, cu, du
   au(:) = 0d0; bu(:) = 0d0; cu(:) = 0d0; du(:) = 0d0

!  set up matrix
   do i=2,N-1
      c     = 2.0*dt*nuY(i)  /(h(i)+h(i+1))/h(i)
      a     = 2.0*dt*nuY(i-1)/(h(i)+h(i-1))/h(i)
!      l     =     dt*Lsour(i)

      cu(i) =-cnpar*c
      au(i) =-cnpar*a
     ! bu(i) = 1d0 + cnpar*(a + c) - l
      bu(i) = 1d0 + cnpar*(a + c)

     ! du(i) = (1d0 - (1d0-cnpar)*(a + c))*Yin(i)                  &
     !       + (1d0 - cnpar)*( a*Yin(i-1) + c*Yin(i+1) ) + dt*Qsour(i)

      du(i) = (1d0 - (1d0-cnpar)*(a + c))*Yin(i)                  &
            + (1d0 - cnpar)*( a*Yin(i-1) + c*Yin(i+1) )

   enddo

!   set up upper boundary condition
!   select case(Bcup)
!   case(Neumann)
      a     = 2.0*dt*nuY(N-1)/(h(N)+h(N-1))/h(N)
!      l     = dt*Lsour(N)

      au(N) =-cnpar*a
!      if (posconc .eq. 1 .and. Yup.lt.0d0) then ! Patankar (1980) trick
!         bu(N) =  1d0 - au(N) - l  - dt*Yup/Yin(N)/h(N)
!         du(N) = Yin(N) + dt*Qsour(N)   &
!               + (1d0 - cnpar)*a*(Yin(N-1)-Yin(N))
!      else
!       bu(N) =  1d0 - au(N) - l
       bu(N) = 1d0 - au(N)
!       du(N) = Yin(N) + dt*(Qsour(N)+Yup/h(N))   &
!             + (1d0 - cnpar)*a*(Yin(N-1)-Yin(N))
       du(N) = Yin(N) + (1d0 - cnpar)*a*(Yin(N-1)-Yin(N))
!      end if
!   case(Dirichlet)
!      au(N) = 0d0
!      bu(N) = 1d0
!      du(N) = Yup
!   case default
!      write(6,*) 'Fatal error: invalid boundary condition type for upper boundary'
!      stop  'diff_center.F90'
!   end select

!   set up lower boundary condition
!   select case(Bcdw)
!   case(Neumann)
      c     = 2.0*dt*nuY(1)/(h(1)+h(2))/h(1)
!      l     = dt*Lsour(1)

      cu(1) =-cnpar*c
     ! if (posconc.eq.1 .and. Ydw.lt.0d0) then ! Patankar (1980) trick
     !    bu(1) = 1d0 - cu(1) - l - dt*Ydw/Yin(1)/h(1)
     !    du(1) = Yin(1) + dt*(Qsour(1))   &
     !          + (1d0 - cnpar)*c*(Yin(2)-Yin(1))
     ! else
     !    bu(1) = 1d0 - cu(1) - l
     !    du(1) = Yin(1) + dt*(Qsour(1)+Ydw/h(1))   &
     !          + (1d0 - cnpar)*c*(Yin(2)-Yin(1))

         bu(1) = 1d0 - cu(1)
         du(1) = Yin(1) + (1d0 - cnpar)*c*(Yin(2)-Yin(1))

     ! endif
!   case(Dirichlet)
!      cu(1) = 0d0
!      bu(1) = 1d0
!      du(1) = Ydw
!   case default
!      write(6,*) 'Fatal error: invalid boundary condition type for lower boundary'
!      stop  'diff_center.F90'
!   end select


!  relaxation to observed value
!   if (minval(Taur).lt.1.E10) then
!      do i=1,N
!         bu(i)=bu(i)+dt/Taur(i)
!         du(i)=du(i)+dt/Taur(i)*Yobs(i)
!      end do
!   end if


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

!   write(6,*) 'fi = ',fi
!   write(6,*) 'value(fi) = ',value(fi)

   do i=fi+1,lt
      value(i)=qu(i)-ru(i)*value(i-1)
   end do

!   write(6,*) 'value(lt) = ',value(lt)

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
End MODULE FlexEFT1D
