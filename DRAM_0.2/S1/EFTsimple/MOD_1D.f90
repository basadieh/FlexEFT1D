MODULE MOD_1D
use BIO_MOD
implicit none
public
! Grid parameters
real, private, parameter :: hmax   = 2.5d2   ! Total water depth
real, private, parameter :: thetaS = 2d0    ! surface stretching parameter
real, private, parameter :: dtsec  = 3d2     ! time step in seconds
real, private, parameter ::d_per_s = 864d2   ! how many seconds in one day
real, private, parameter :: Yup    = 0d0
real, private, parameter :: Ydw    = 0d0 
                      !how many seconds of one year
integer, private, parameter :: y_per_s = INT(d_per_s*360), &

!  Number of vertical points of NO3,Temp, and Aks profile
                      N_NO3   = 37,               &  
                      N_Temp  = 57,               &
                      N_w     = 40,               &
                      N_par   = 1 ,               &
                      nsave   = INT(d_per_s)/INT(dtsec) ! Timesteps to save
integer, private            :: N_Aks                    ! Station dependent
!  Indices for external forcing variables
integer, private, parameter :: etemp      = 1
integer, private, parameter :: eNO3       = 2
integer, private, parameter :: eAks       = 3
integer, private, parameter :: ew         = 4
integer, private, parameter :: ePAR       = 5
integer, private, parameter :: TNFo       = ePAR  ! TOtal number of forcings

! Total of observation times in forcing data
! Use Taketo's Aks data for K2 and my own ROMS data for S1 
character(LEN=5), private, parameter :: LabelForc(TNFo) = (/'temp','NO3','Aks','wROMS','par'/)
integer,          private, parameter :: NFobs(TNFo)     = (/    12,   12, 360,     12,  12 /)

! Forcing data time indices 
real, private :: obs_time_temp(NFobs(etemp))
real, private :: obs_time_NO3( NFobs(eNO3))
real, private :: obs_time_Aks( NFobs(eAks))
real, private :: obs_time_w(   NFobs(ew  ))
real, private :: obs_time_par( NFobs(ePAR))

! Forcing data
real, private :: obs_NO3( N_NO3 ,  1+NFobs(eNO3 ))
real, private :: obs_Temp(N_Temp,  1+NFobs(etemp))
real, private :: obs_PAR( N_par ,  1+NFobs(ePAR ))
real, private :: obs_w(   N_w   ,  1+NFobs(ew   ))
real, private, allocatable :: obs_Aks(:,:) ! (N_Aks ,  1+NFobs(eAks ))
 ! Vertically interpolated temperature and Aks at each obs timing
real, private :: VTemp(nlev, NFobs(etemp))
real, private :: VAks(0:nlev,NFobs(eAks ))
real, private :: Vw(  0:nlev,NFobs(ew ))

logical, public  :: savefile

! New calculated state variables
real, private, allocatable    :: ww(:,:)  

!!$$-----------------------------------------------------------------
!$ The declaration of the data part:
! length of string for labels 
integer, parameter   :: LabelLen = 7 

! NDTYPE is the number of types of data (or observations), 
integer              :: NDTYPE
character(LabelLen), allocatable :: DataLabel(:)

! Number of observations (data points) of each data type
integer, allocatable :: NDPTS(:)
integer              :: TNobs

! Indeces for data type:
integer, parameter   :: itNO3 = 1, itCHL = 2, itNPP = 3, itP10 = 4 
integer, parameter   :: itP03 = 5, itP01 = 6, itP_1 = 7

! The number of days for model run, needs to change at different stages during MCMC
integer              :: NDays

! The data arrays are read from csv files. 
! Dimensions of the input data
integer, allocatable :: ncol(:)
integer, allocatable :: nrow(:)

! Declare the dimensions of the data:
real,    allocatable ::  TINData(:,:)
real,    allocatable ::  CHLData(:,:)
real,    allocatable ::  NPPData(:,:)
real,    allocatable :: SizeData(:,:)
real,    allocatable ::  OBSData(:,:)

! DOY and Depth for Observational data assembled as a single matrix 
real,    allocatable ::  OBS_DOY(:), OBS_Depth(:)
! Data label for OBS data
character(LabelLen), allocatable :: OBS_Label(:)

! Initial profile of NO3:
real, private               :: NO3(nlev,1)

CONTAINS
!-----------------------------------------------------------------------
! Call other subroutines to read and set-up data, as needed 
! This subroutine is called before the real model run
subroutine Setup_OBSdata
implicit none
character(LEN=20)    :: TIN_OBS_file, CHL_OBS_file, NPP_OBS_file, SIZE_OBS_file
real,    allocatable :: DOY(:), Depth(:)
integer              :: k, i,oi

  write(6,*) 'Stn name: ',Stn
  Select case(Model_ID)
  case(Geidersimple, EFTsimple)

     NDTYPE = 3  !TIN, CHL, PP
     allocate(DataLabel(NDTYPE),  STAT = AllocateStatus)
     IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

     DataLabel(itNO3) = 'TIN'
     DataLabel(itCHL) = 'CHL'
     DataLabel(itNPP) = 'NPP'

     allocate(NDPTS(NDTYPE),  STAT = AllocateStatus)
     IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

! Assign the data dimension, must be consistent with external file:
     if (trim(Stn) .eq. 'S1') then
         N_Aks        = 40
         NDPTS(itNO3) = 902
         NDPTS(itCHL) = 1993
         NDPTS(itNPP) = 128

     else if (trim(Stn) .eq. 'K2') then

         N_Aks        = 25
         NDPTS(itNO3) = 974
         NDPTS(itCHL) = 1253
         NDPTS(itNPP) = 112

     else if (trim(Stn) .eq. 'HOT') then
         N_Aks        = 40
         NDPTS(itNO3) = 3910
         NDPTS(itCHL) = 8180
         NDPTS(itNPP) = 1659

     else
         write(6,*) 'Station number incorrect! Stop!'
         stop
     endif

     TNobs  = sum(NDPTS(:))

     allocate(nrow(3), STAT = AllocateStatus)
     IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

     allocate(ncol(3), STAT = AllocateStatus)
     IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

     do i = 1,3
        nrow(i)=NDPTS(i)
        ncol(i)=3
     enddo

  case(EFTdiscrete,EFTcont)

     NDTYPE = 7  !TIN, CHL, PP, CHL>10, CHL3-10,CHL1-3,CHL<1
     allocate(DataLabel(NDTYPE),  STAT = AllocateStatus)
     IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

     DataLabel(itNO3) = 'TIN'
     DataLabel(itCHL) = 'CHL'
     DataLabel(itNPP) = 'NPP'
     DataLabel(itP10) = 'P10'
     DataLabel(itP03) = 'P03'
     DataLabel(itP01) = 'P01'
     DataLabel(itP_1) = 'P_1'

     allocate(NDPTS(NDTYPE),  STAT = AllocateStatus)
     IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

     allocate(nrow(4), STAT = AllocateStatus)
     IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

     allocate(ncol(4), STAT = AllocateStatus)
     IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

     if (trim(Stn) .eq. 'S1') then

         N_Aks        = 40
         NDPTS(itNO3) = 902
         NDPTS(itCHL) = 1993
         NDPTS(itNPP) = 128
         NDPTS(itP10) = 169

     else if (trim(Stn) .eq. 'K2') then

         N_Aks        = 25
         NDPTS(itNO3) = 974
         NDPTS(itCHL) = 1253
         NDPTS(itNPP) = 112
         NDPTS(itP10) = 143

     else if (trim(Stn) .eq. 'HOT') then
         write(6,*) 'No size data available for HOT! Quit!'
         stop
     else
         write(6,*) 'Station number incorrect! Stop!'
         stop
     endif

     NDPTS(itP03) = NDPTS(itP10)
     NDPTS(itP01) = NDPTS(itP10)
     NDPTS(itP_1) = NDPTS(itP10)
     TNobs        = sum(NDPTS(:))

     do i = 1,4
        nrow(i)=NDPTS(i)
        if (i < 4) then
           ncol(i)=3
        else
           ncol(i)=6
        endif
     enddo
  case default
     write(6,*) 'Model option incorrect! Quit!'
     stop

  End select

  write(6,'(A30,1x,I5)') 'Total number of observations =',TNobs
 
  allocate( TINData(nrow(1),ncol(1)), STAT = AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

  allocate( CHLData(nrow(2),ncol(2)), STAT = AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

  allocate( NPPData(nrow(3),ncol(3)), STAT = AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

  allocate( OBSData(TNobs  ,3      ), STAT = AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  ! Initialize observational data:
  TINData (:,:) = 0d0
  CHLData (:,:) = 0d0
  NPPData (:,:) = 0d0
  OBSData (:,:) = 0d0

  ! Read observational data:
  TIN_OBS_file = trim(Stn)//'_TIN.dat' 
  CHL_OBS_file = trim(Stn)//'_CHL.dat'
  NPP_OBS_file = trim(Stn)//'_NPP.dat'
  
  call Readcsv( TIN_OBS_file,nrow(1),ncol(1), TINData)
  call Readcsv( CHL_OBS_file,nrow(2),ncol(2), CHLData)
  call Readcsv( NPP_OBS_file,nrow(3),ncol(3), NPPData)

  if (Model_ID == EFTdiscrete .or. Model_ID == EFTcont) then
     allocate(SizeData(nrow(4),ncol(4)), STAT = AllocateStatus)
     IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

     SizeData(:,:) = 0d0
     SIZE_OBS_file = trim(Stn)//'_size.dat'
     call Readcsv(SIZE_OBS_file,nrow(4),ncol(4),SizeData)
  endif
 

  allocate(OBS_DOY(  TNobs))
  allocate(OBS_Depth(TNobs))
  allocate(OBS_Label(TNobs))

  OBS_DOY(:)   = 0d0
  OBS_Depth(:) = 0d0
  OBS_Label(:) = ' NA    '

  ! Give DOY and Depth for all the obs. data
  k = 0
  do i = 1, NDTYPE
    allocate(  DOY(NDPTS(i)))
    allocate(Depth(NDPTS(i)))

    selectcase(i)
    
    case(itNO3)   ! Nitrate
      DOY   = TINData(:,1)
      Depth = TINData(:,2)
    case(itCHL)   ! CHL
      DOY   = CHLData(:,1)
      Depth = CHLData(:,2)
    case(itNPP)   ! NPP
      DOY   = NPPData(:,1)
      Depth = NPPData(:,2)

    case(itP10,itP03,itP01,itP_1)   ! Size-fractionated Chl
      DOY   = SizeData(:,1)
      Depth = SizeData(:,2)

    case default
      print *, 'Errors in selecting data types! Quit!'
      stop
    endselect

    do oi = 1, NDPTS(i)
         OBS_DOY(k+oi) =   DOY(oi)
       OBS_Depth(k+oi) = Depth(oi)
       OBS_Label(k+oi) = DataLabel(i)
    enddo
    k = k + NDPTS(i)
    deallocate(DOY)
    deallocate(Depth)
  enddo
  write(6,*) 'Observation data set up!'
end subroutine Setup_OBSdata
!========================================================
subroutine Model_setup
! Called at the beginning of the program, only once
implicit none
integer            :: i
! Input files
character(LEN=20)  :: forcfile(TNFo), forcfiletime(TNFo)

!the fraction of a time step in one day
dtdays = dtsec/d_per_s
write(6,'(A13,1x,F12.3,A12)') 'Timestepping: ',dtdays,'of one day.'
! Setup grid
call setup_grid

! Forcing files:
do i = 1, TNFo
   forcfile(i)     = trim(Stn)//'_'//trim(LabelForc(i))//'.dat'
   forcfiletime(i) = trim(Stn)//'_'//trim(LabelForc(i))//'_time.dat'
enddo

! Use Taketo's Aks data only for K2:
if (trim(Stn) .eq. 'K2') then
   forcfile(eAks)     = trim(Stn)//'_'//trim(LabelForc(eAks))//'T.dat'
   forcfiletime(eAks) = trim(Stn)//'_'//trim(LabelForc(eAks))//'_timeT.dat'
endif

! Read NO3 data:
call Readcsv(forcfile(eNO3), N_NO3, size(obs_NO3,2), obs_NO3) 

! Interpolate initial NO3:
call gridinterpol(N_NO3,1,obs_NO3(:,1),obs_NO3(:,2),   &
                  nlev, Z_r, NO3(:,1)) 

! Calculate obs. time indices in seconds:
! Read obs. time file:
call Readcsv(forcfiletime(ePAR), 1,NFobs(ePAR), obs_time_par) 
call Readcsv(forcfiletime(eNO3), 1,NFobs(eNO3), obs_time_NO3) 
call Readcsv(forcfiletime(etemp),1,NFobs(etemp),obs_time_temp) 
call Readcsv(forcfiletime(eAks), 1,NFobs(eAks), obs_time_Aks) 
call Readcsv(forcfiletime(ew  ), 1,NFobs(ew)  , obs_time_w) 

! Convert obs_time to the unit of seconds:
obs_time_par  = obs_time_par *3d1*dble(d_per_s)
obs_time_w    = obs_time_w   *3d1*dble(d_per_s)
obs_time_Aks  = obs_time_Aks *3d1*dble(d_per_s)
obs_time_temp = obs_time_temp*3d1*dble(d_per_s)
obs_time_NO3  = obs_time_NO3 *3d1*dble(d_per_s)

! Read external w data:
!call Readcsv(forcfile(ew), size(obs_w,1), size(obs_w,2), obs_w) 

! Interpolate external w data:

!subroutine gridinterpol(N,cols,obs_z,obs_prof,nlev_,model_z,model_prof)
!call gridinterpol(N_w,NFobs(ew),obs_w(:,1),       &
!      obs_w(:,2:(NFobs(ew)+1)),                   &
!      nlev+1, Z_w, Vw(:,:)) 

! Read external PAR data:
call Readcsv(forcfile(ePAR), size(obs_PAR,1),size(obs_PAR,2),obs_PAR) 

! Read external Temp data:
call Readcsv(forcfile(etemp),size(obs_Temp,1),  &
             size(obs_Temp,2),obs_Temp) 

! Interpolate external Temp data:
call gridinterpol(N_Temp,NFobs(etemp),obs_Temp(:,1),      &
      obs_Temp(:,2:(NFobs(etemp)+1)),                     &
      nlev, Z_r, VTemp(:,:)) 

! Read external Aks data:
allocate(obs_Aks(N_Aks, 1+NFobs(eAks)),  STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
call Readcsv(forcfile(eAks), N_Aks, size(obs_Aks,2), obs_Aks) 

! Interpolate external Aks data:
call gridinterpol(N_Aks,NFobs(eAks),obs_Aks(:,1), &
      obs_Aks(:,2:(NFobs(eAks)+1)),               &
      nlev+1, Z_w, VAks(:,:)) 

savefile = .FALSE.
write(6,*) 'Model setup finished!'
End subroutine Model_setup
!========================================================
! This subroutine does the main work in the 1D model
! and give the output to match with the observational data
! Must be called after initialization
SUBROUTINE Timestep(TINout, CHLout, NPPout, Sizeout)
implicit none
real,    parameter  :: cnpar = 0.6

! Local scratch variables
integer  :: it,i,k, j, Nstep, current_day, current_DOY,DOY
real     :: cff,current_sec
real     :: I_0(1),Aks(0:nlev),w(0:nlev)
real     :: a(nlev,1),depth(1)
real     :: CHL_(nlev),Vars1(nlev),Vars2(nlev),ww_(0:nlev)
real     :: VPAR(NFobs(ePAR))

integer, parameter :: mode0 = 0
integer, parameter :: mode1 = 1
real,    parameter :: Aks_th= 1D-3 ! Aks threshold for calculating average PAR in MLD 

character(LEN=20)  :: outfile
! The model output (final year) to match with observational data: 
real,  intent(out)           ::  TINout(size( TINData,1),1)
real,  intent(out)           ::  CHLout(size( CHLData,1),1)
real,  intent(out)           ::  NPPout(size( NPPData,1),1)
real,  intent(out), optional :: Sizeout(size(SizeData,1),4)
logical                      :: MLD_found
! Initialize initial NO3:
Vars(iNO3,:) = NO3(:,1)

! Initialize other variables:
do k = 1,nlev
   do i = 1,NPHY
      Vars(iPHY(i), k) = 0.1/float(NPHY)
   enddo
   Vars(iZOO,k) = 0.1
   Vars(iDET,k) = 0.1

   if (NVAR > iDET) then
      do i = (iDET+1), NVAR
         Vars(i,k) = 1D-1
      enddo
   endif
enddo

! Initialize output data:
Varout(:,:) = 1D-3
do i = 1,NVAR
   do k = 1,nlev
      Varout(i,k)= Vars(i,k)
   enddo
enddo

 TINout(:,:) = 0d0
 CHLout(:,:) = 0d0
 NPPout(:,:) = 0d0

if (Model_ID == EFTdiscrete .or. Model_ID == EFTcont) then
   Sizeout(:,:) = 0d0
endif

! Sinking rate
allocate(ww(0:nlev,NVsinkterms),  STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

! Create output files:
if (savefile) then
   do i = 1, (Nout+ow)
      outfile = trim(Labelout(i))//'.out'

      ! Create data out file
      open (9+i, file = outfile, status = 'replace')
    
      if ( (i .eq. oAks) .or. (i .eq. ow) ) then
         write(9+i, 3) 'Timestep','Days',Z_w
      else
         write(9+i, 1) 'Timestep','Days',Z_r
      endif
   enddo
endif

1 format(A10,1x,A7, <nlev  >(2x,F12.5))
3 format(A10,1x,A7, <nlev+1>(2x,F12.5))

! Number of time steps
Nstep = NDays*INT(d_per_s)/INT(dtsec) 

! 'START TIME STEPPING'
DO it = 1, Nstep+1

! Calculate current timing (zero is starting time):
   current_sec = float(it-1)*dtsec

! For each time step, read in external environmental data
! Interpolate Temp data throughout the water column:
!subroutine time_interp(time, N, nrow, obs_time, obs_data, mod_data)

   call time_interp(int(current_sec),size(obs_time_temp,1),&
                    nlev,obs_time_temp,VTemp,Temp)

  ! Interpolate temporal PAR data: 
   VPAR = obs_PAR(1,2:size(obs_PAR,2))
   call time_interp(int(current_sec),size(obs_time_par,1),1,&
        obs_time_par,VPAR,I_0)

  ! Convert the unit of par to W m-2:
   cff = I_0(1)/4d-1

  ! Calculate light field:
   CHL_(:) = Varout(oCHLt,:)
   call Calculate_PAR(cff, nlev, Hz, CHL_, PAR)
     
  ! Interpolate Aks data throughout the water column:
   call time_interp(int(current_sec), size(obs_time_Aks,1), nlev+1,&
        obs_time_Aks, VAks, Aks)

  ! Calculate the vertical grid index (N_MLD) at the bottom of MLD:
   MLD_found = .FALSE.
   do k = nlev, 1, -1
      if (Aks(k) .gt. Aks_th .and. Aks(k-1) .le. Aks_th) then
         N_MLD     = k
         MLD_found = .TRUE.
         exit
      endif
   enddo
   if (.not. MLD_found) N_MLD=1  ! Mixed throughout the whole water column

  ! Calculate average PAR within the surface mixed layer (from nlev to N_MLD):
   PARavg=0d0
   do k=nlev,N_MLD,-1
      PARavg=PARavg+PAR(k)*Hz(k)
   enddo 
   PARavg = PARavg/abs(Z_w(N_MLD-1))

   !write(6,*) 'N_MLD = ',N_MLD
   !write(6,*) 'MLD   = ',Z_w(N_MLD-1)
   !write(6,*) 'PARavg = ',PARavg
   !stop
  ! Interpolate w data throughout the water column:
  ! call time_interp(int(current_sec), size(obs_time_w,1), nlev+1,&
  !      obs_time_w, Vw, w)
  ! w(0   ) = 0d0
  ! w(nlev) = 0d0

  IF (mod(it, nsave) .EQ. 1) THEN

    ! Check whether the values are valid:
    do j = 1,NVAR
       do k=1,nlev
          if( (Vars(j,k) .ne. Vars(j,k)) .OR. (Vars(j,k) .le. 0.) ) then
              write(6,*) 'j = ',j
              write(6,*) 'k = ',k
              write(6,*) 'At time step ',it
              write(6,*) 'The variable ',trim(Labelout(j+ow)), &
                        ' is invalid at depth ',Z_r(k)
              write(6,*) 'Vars(j,k) =',  Vars(j,k)
              stop
          endif
       enddo 
    enddo

   ! Calculate model time in days:
     current_day = int(current_sec/d_per_s)

  ! Calculate DATE OF the YEAR (DOY)
     current_DOY = mod(current_day, 360)

    ! Save data to output files:
    if (savefile) then

    ! Save Temp data into the output file:
       write(9+oTemp, 2) it, current_day, Temp

    ! Save PAR data into the output file:
       write(9+oPAR, 2) it, current_day, PAR

    ! Save Aks data into the output file:
       write(9+oAks, 4) it, current_day, Aks

    ! Save w data into the output file:
    !   write(9+ow  , 4) it, current_day, w

    ! Save state variables and diagnostics:
       do i=(ow+1),(Nout+ow)
          Vars1(:) = Varout(i-ow,:)
          write(9+i, 2) it, current_day, Vars1
       enddo

    endif  !==> End of saving results
  
!-------------------------------------------------------------------------
  ! Calculate model outputs (final year) to match with obs. data
 If ((NDays-current_day) .le. 360) Then

    Do i=1, sum(nrow(:))
       DOY      =  INT(min(OBS_DOY(i),360.0))
       depth(1) = -abs(OBS_Depth(i))  !Consistent with ROMS convention

       IF (DOY .eq. current_DOY) then

           if (i .le. nrow(1)) then
           ! Calculate TIN output:
           
            a(:,1) = Vars(iNO3,:)
            call gridinterpol(nlev,1,Z_r(:),a(:,1),1,depth(1),&
                 TINout(i,1))       
           
        
           elseif (i .le. (nrow(1)+nrow(2)) ) then

           ! Calculate CHL output:
            a(:,1) = Varout(oCHLt,:)
         
            call gridinterpol(nlev,1,Z_r(:),a(:,1),1,depth(1),&
                 CHLout((i-nrow(1)),1))   

           elseif (i .le. (nrow(1)+nrow(2)+nrow(3)) ) then

           ! Calculate NPP output:
            a(:,1) = Varout(oPPt,:)   ! Carbon based PP
            call gridinterpol(nlev,1,Z_r(:),a(:,1),1,depth(1),&
                 NPPout((i-nrow(1)-nrow(2)),1))  

           else
             if (Model_ID == EFTdiscrete .or. Model_ID == EFTcont) then
           !    Calculate Size-fractionated Chl output:
                do j = 1, 4
                   a(:,1) = Varout(oCHLs(j),:) 

                   call gridinterpol(nlev,1,Z_r(:),a(:,1),1,depth(1), &
                        Sizeout((i-nrow(1)-nrow(2)-nrow(3)),j))
                enddo
             endif
           endif
       Endif !==> End of matching DOY
    Enddo
 Endif !==> End of checking if final year
!---------------------------------------------------------------------
 ENDIF  !==> END of daily work
2 format(I10,1x,I7, <nlev  >(2x,1pe12.3))
4 format(I10,1x,I7, <nlev+1>(2x,1pe12.3))
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
    call FlexEFT_simple
  case(Geidersimple)
    call Geider_simple
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

   Vars1(:)     = Vars(j,:)
   call diff_center(nlev,dtsec,cnpar,1,Hz,Aks,Vars1,Vars2)

   ! Save diffusion fluxes (normalized to per day)
   Varout(oD_NO3+j-1,:) = (Vars2(:) - Vars1(:))/dtdays

   ! Update the state variables:
   Vars(j,:)=Vars2(:)
  
enddo

! Sinking:
! Initialize sinking rate (UNIT: m/s !):
ww(0,   :) = 0d0
ww(nlev,:) = 0d0

do k = 1,nlev-1
   do i=1,NPHY
      !Phytoplankton sinking rate  
      ww(k,i) = -Varout(ow_p(i),k)/dble(d_per_s)
   enddo
 !Detritus sinking rate (convert to UNIT: m/s)
   ww(k,NPHY+1) = -Params(iwDET)/dble(d_per_s) 
enddo

do j = 1,NVsinkterms
   ww_(:)   = ww(:,j)
   Vars2(:) = Vars(Windex(j),:)
   call adv_center(nlev,dtsec,Hz,Hz,ww_(:),1,1,Yup,Ydw,6,mode1,Vars2(:))
   Vars(Windex(j),:) = Vars2(:)
enddo

! Vertical advection (due to w, unit: m/s):
! Now do not consider vertical w
!do j = 1,NVAR
!   call adv_center(nlev,dtsec,Hz,Hz,w(:),  &
!                   1,1,Yup,Ydw,6,mode0,Vars(j,:))
!enddo

!write(6,*) 'After advection, Vars(iCHL,1) = ',Vars(iCHL,1)

ENDDO    ! ==> End of time stepping

! Close files:
if (savefile) then
   do i = 1, (Nout+ow)
      close (unit=9+i)
   enddo
endif
deallocate(ww)
END subroutine Timestep
!=====================================================
subroutine End_model
  implicit none

!  ! Release memory:
  deallocate(Vars)
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

END subroutine End_model
!========================================================
subroutine setup_grid
implicit none

! Local scratch variable:
integer                :: i
real                   :: sc_r, C_sig

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
end subroutine setup_grid
!===========
subroutine Calculate_PAR(I_0, nlev_, Hz, Chl, PAR)
  !top level is nlev_, bottom layer is 1, following ROMS convention
  implicit none
  real   , intent(in) :: I_0
  integer, intent(in) :: nlev_    ! Total number of vertical layers
  real   , intent(in) :: Hz(nlev_), Chl(nlev_)
  real   , intent(out):: PAR(nlev_)
  integer             :: i
  real                :: par0, attn    ! Scratch variable
  real   , parameter  :: kw = 0.04, kc = 0.025

  par0 = I_0   !Light at the grid surface
  do i = nlev_,1,-1
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
   real,     intent(in)                 :: dt
!  layer thickness (m)
   real,     intent(in)                 :: h(N)

!  old layer thickness (m)
   real,     intent(in)                 :: ho(N)

!  vertical advection speed (m/s)
   real,     intent(in)                 :: ww(0:N)

!  type of upper BC
   integer,  intent(in)                :: Bcup

!  type of lower BC
   integer,  intent(in)                :: Bcdw

!  value of upper BC
   real, intent(in)                 :: Yup

!  value of lower BC
   real,     intent(in)            :: Ydw

!  type of advection scheme
   integer,  intent(in)            :: method

!  advection mode (0: non-conservative, 1: conservative)
   integer,  intent(in)            :: mode
!
! !INPUT/OUTPUT PARAMETERS:
   real,  intent(inout)            :: Y(N)
!
! !DEFINED PARAMETERS:
   real,     parameter             :: one6th=1d0/6d0
   integer,  parameter             :: itmax=100

!  type of advection scheme
   integer,  parameter             :: UPSTREAM       = 1
   integer,  parameter             :: P1             = 2
   integer,  parameter             :: P2             = 3
   integer,  parameter             :: Superbee       = 4
   integer,  parameter             :: MUSCL          = 5
   integer,  parameter             :: P2_PDM         = 6

!  boundary condition type
!  for advection schemes
   integer,  parameter             :: flux           = 1
   integer,  parameter             :: value          = 2
   integer,  parameter             :: oneSided       = 3
   integer,  parameter             :: zeroDivergence = 4

! !LOCAL VARIABLES:
   integer                         :: i,k,it
   real                            :: x,r,Phi,limit
   real                            :: Yu,Yc,Yd
   real                            :: c,cmax
   real                            :: cu(0:N)
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

!   limit the flux according to different suggestions (Pietrzak 1998)
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
            Y(k)=Y(k)-1d0/float(it)*dt*((cu(k)-cu(k-1))/        &
                 h(k)-Y(k)*(ww(k)-ww(k-1))/h(k))
         enddo
      else                ! conservative, for PHY and DET sinking
         do k=1,N
            Y(k)=Y(k)-1d0/float(it)*dt*((cu(k)-cu(k-1))/h(k))
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
End MODULE MOD_1D
