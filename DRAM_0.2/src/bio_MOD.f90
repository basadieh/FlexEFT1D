Module BIO_MOD
implicit none
private

! !PUBLIC MEMBER FUNCTIONS:
public choose_model,FLEXEFT_DISC,FLEXEFT_CONT,Cal_TN
public Geider_simple, FlexEFT_simple
                                ! Number of vertical layers
integer, public, parameter :: nlev        = 30,       &  
                              ! Options of biological models
                              EFTdiscrete = 1,        &
                              EFTcont     = 2,        &
                              Geiderdisc  = 3,        &
                              Geidercont  = 4,        &
                              EFTsimple   = 5,        &
                              Geidersimple= 6

! Parameters for phytoplankton size fractional Chl
real, public, parameter :: pi=3.1415926535897932384633D0
real, parameter :: PMU_min = log(1d1*pi/6d0*0.5**3)
real, parameter :: PMU_max = log(1d1*pi/6d0*6d1**3)
real, parameter :: PMU_1   = log(1d1*pi/6d0)
real, parameter :: PMU_3   = log(1d1*pi/6d0*3d0**3)
real, parameter :: PMU_10  = log(1d1*pi/6d0*1d1**3)

integer, public  :: AllocateStatus
integer, public  :: N_MLD   ! vertical grid index at the bottom of MLD


real, public     :: Temp(nlev), PAR(nlev), dtdays, Ntot, PARavg 
real, public     :: Z_r(1:nlev), Z_w(0:nlev), Hz(nlev)  ! Grid variables

integer, public  :: NVAR, Nout, iZOO, iDET,  iPMU, iVAR, iCHL
integer, public  :: NVsinkterms,NPHY, NPar
integer, public  :: oZOO, oDET, oFER, oZ2N, oD2N, oPHYt,oCHLt,oPPt
integer, public  :: oPMU, oVAR, odmudl,od2mu,od2gdl  
integer, public  :: oD_NO3,oD_ZOO,oD_DET,oD_PMU,oD_VAR,oD_CHL
integer, public  :: oCHLs(4)   ! Four size fractions of CHL

! Indices for parameters used in DRAM
integer, public  :: imu0,iaI0,igmax,iKN, ithetm
!integer, public  :: iPHYini, iZOOini, iDETini
integer, public  :: ialphaI,iA0N,ialphaA,ialphaG,ialphaK
integer, public  :: ikp,iQ0N,ialphaQ,iPenfac,iLref,iwDET,irdN,imz
integer, public  :: ithetamin,iQNmin
integer, public, allocatable :: iPHY(:),oPHY(:),oTheta(:),oQN(:),    &
                                omuNet(:),ow_p(:),oD_PHY(:),         &
                                oGraz(:), oSI(:), oLno3(:), oTheHat(:) 
                                             
integer, public, allocatable :: Windex(:)
real,    public, allocatable :: Vars(:,:),Varout(:,:), params(:)

character(LEN=6), public, allocatable :: ParamLabel(:)

! Fixed Model parameters:
real, parameter :: PMUmax =2D1, VARmax=1D2
real, parameter :: w_p0   =0d0, alphaW=0d0, Wmax=5d1
real, parameter :: Ep     =0.4, Ez    =0.6
real, parameter :: Femin  =0.02,K0Fe  =0.8, alphaFe=0.14
real, parameter :: zetaChl=0.8, zetaN =0.6, RMchl0 =0.1
real, parameter :: GGE    =0.3, unass =0.24
real, parameter :: alphamu=0.0, betamu=0.0
real, parameter :: thetm  =0.65
real, parameter :: mu0    =5d0, V0N   =5d0
real, parameter :: alphaV =0d0

! Size and weights for discrete models
real, allocatable          :: PMU_(:)
real, allocatable          :: wtCHL(:,:)  ! weight for each size class

! Indices for state variables
integer, public, parameter :: iNO3=1,oTemp=1,oPAR=2,oAks=3,oNO3=1,ow=4
integer, public            :: nutrient_uptake=1, grazing_formulation=3   

! Output indices:
character(LEN=10), public, allocatable  :: Labelout(:)  

integer, public,   parameter :: namlst=8
character(LEN=3),  public    :: Stn

! Model options:
integer,           public    :: Model_ID = EFTdiscrete

! Local variables:
real :: INGES,gbar,EGES,Zmort,RES
real :: PHY,NO3,ZOO,DET,DET1
CONTAINS
!========================================================
subroutine choose_model
implicit none
integer             :: i
namelist /station/    Stn, Model_ID, nutrient_uptake, grazing_formulation

!  open the namelist file and read station name.
open(namlst,file='Station.nml',status='old',action='read')
read(namlst,nml=station)
close(namlst)

SELECT CASE(model_ID)
  case(Geidersimple)
    write(6,*) 'Geider simple model selected!'
    NPHY = 1
    allocate(iPHY(NPHY))

    do i=1,NPHY
       iPHY(i)=i+iNO3
    enddo 

    iZOO=iPHY(NPHY)+1
    iDET=iZOO+1
    iCHL=iDET+1
    NVAR=iCHL 

    allocate(Vars(NVAR,nlev))

    NVsinkterms = 1+NPHY
    allocate(Windex(NVsinkterms))
    do i=1,NPHY
       Windex(i)=iPHY(i)
    enddo
    Windex(NVsinkterms)=iDET

    ! Output array matrices
    allocate(oPHY(NPHY))
    do i=1,NPHY
       oPHY(i)=i+oNO3
    enddo
    oZOO =oPHY(NPHY)+1
    oDET =oZOO+1
    oCHLt=oDET+1
    ! The above must match with i** indeces   
 
    allocate(omuNet(NPHY))
    allocate(oGraz(NPHY))
    allocate(ow_p(NPHY))
    allocate(oD_PHY(NPHY))

    do i=1,NPHY
       omuNet(i)=oCHLt + i
    enddo

    do i=1,NPHY
       oGraz(i)=omuNet(NPHY)+i
    enddo

    oZ2N=oGraz(NPHY)+1
    oD2N=oZ2N+1
    do i=1,NPHY
       ow_p(i)=oD2N+i
    enddo


    oPPt  =ow_p(NPHY)+1
    oD_NO3=oPPt +1

    do i=1,NPHY
       oD_PHY(i)=oD_NO3+i
    enddo

    oD_ZOO=oD_PHY(NPHY)+1
    oD_DET=oD_ZOO+1
    oD_CHL=oD_DET+1
    Nout  =oD_CHL

    allocate(Varout(  Nout,nlev))
    allocate(Labelout(Nout+ ow ))

    Labelout(oTemp  )='Temp '
    Labelout(oPAR   )='PAR  '
    Labelout(oAks   )='Aks  '
    Labelout(ow     )='w    '
    Labelout(oNO3+ow)='NO3  '

    do i=1,NPHY
       write(Labelout(oPHY(i)+ow), "(A3,I2)") 'PHY',i
    enddo

    Labelout(oZOO +ow)='ZOO  '
    Labelout(oDET +ow)='DET  '
    Labelout(oCHLt+ow)='CHL_T'
    do i=1,NPHY
       write(Labelout(omuNet(i) + ow), "(A3,I2)") 'muN',i
       write(Labelout(oGraz(i)  + ow), "(A3,I2)") 'Gra',i
       write(Labelout(ow_p(i)   + ow), "(A3,I2)") 'w_p',i
       write(Labelout(oD_PHY(i) + ow), "(A3,I2)") 'D_P',i
    enddo
  
    Labelout(oZ2N  + ow)='Z2N  '
    Labelout(oD2N  + ow)='D2N  '
    Labelout(oPPt  + ow)='NPP_T'
    Labelout(oD_NO3+ ow)='D_NO3'
    Labelout(oD_ZOO+ ow)='D_ZOO'
    Labelout(oD_DET+ ow)='D_DET'
    Labelout(oD_CHL+ ow)='D_CHL'

    do i = 1, Nout+ow
       write(6,*) 'Labelout(',i,') = ',Labelout(i)
    enddo
    ! Initialize parameters
    ! Indices for parameters that will be used in MCMC                 
    imu0    =  1
    iaI0    =  imu0   + 1
    iQ0N    =  iaI0   + 1
    iKN     =  iQ0N   + 1
    igmax   =  iKN    + 1
    ikp     =  igmax  + 1 
    iwDET   =  ikp    + 1
    irdN    =  iwDET  + 1
    imz     =  irdN   + 1
 !   iPHYini =  imz    + 1
 !   iZOOini =  iPHYini+ 1
 !   iDETini =  iZOOini+ 1
 !   NPar    =  iDETini
    NPar    =  imz
    write(6,'(I2,1x,A20)') NPar,'parameters in total to be estimated.'
    allocate(params(NPar))
    allocate(ParamLabel(NPar))

    ParamLabel(imu0  ) = 'mu0hat '
    ParamLabel(iaI0  ) = 'aI0    '
    ParamLabel(iKN   ) = 'KN     '
    ParamLabel(iQ0N  ) = 'Q0N    '
    ParamLabel(iwDET ) = 'wDET   '
    ParamLabel(igmax ) = 'gmax   '
    ParamLabel(ikp   ) = 'kp     '
    ParamLabel(irdn  ) = 'rdn    '
    ParamLabel(imz   ) = 'mz     '
  !  ParamLabel(iPHYini)= 'PHYini '
  !  ParamLabel(iZOOini)= 'ZOOini '
  !  ParamLabel(iDETini)= 'DETini '

    params(imu0   ) = 2d0
    params(iaI0   ) = 0.5
    params(iKN    ) = 1d0
    !params(iPHYini) = 0.1
    !params(iZOOini) = 0.1
    !params(iDETini) = 0.1
    params(igmax  ) = 1d0
    params(ikp    ) = 5d-1
    params(iQ0N   ) = 0.15
    params(iwDET  ) = 1D0
    params(irdn   ) = 5d-2
    params(imz    ) = 5d-2

   CASE(EFTsimple)

    write(6,*) 'FlexEFT simple model (only ONE PHYTO size class) selected!'
    NPHY = 1
    allocate(iPHY(NPHY))

    do i=1,NPHY
       iPHY(i)=i+iNO3
    enddo 

    iZOO=iPHY(NPHY)+1
    iDET=iZOO+1
    NVAR=iDET 

    allocate(Vars(NVAR,nlev))

    NVsinkterms = 1+NPHY
    allocate(Windex(NVsinkterms))
    do i=1,NPHY
       Windex(i)=iPHY(i)
    enddo
    Windex(NVsinkterms)=iDET

    ! Output array matrices
    allocate(oPHY(NPHY))
    do i=1,NPHY
       oPHY(i)=i+oNO3
    enddo
    oZOO =oPHY(NPHY)+1
    oDET =oZOO+1
    oCHLt=oDET+1

    ! The above must match with i** indeces   
 
    allocate(omuNet(NPHY))
    allocate(oGraz(NPHY))
    allocate(oTheta(NPHY))
    allocate(oQN(NPHY))
    allocate(oSI(NPHY))
    allocate(oLno3(NPHY))
    allocate(ow_p(NPHY))
    allocate(oD_PHY(NPHY))

    do i=1,NPHY
       omuNet(i)=oCHLt + i
    enddo

    do i=1,NPHY
       oGraz(i)=omuNet(NPHY)+i
    enddo

    oZ2N=oGraz(NPHY)+1
    oD2N=oZ2N+1
    do i=1,NPHY
       ow_p(i)=oD2N+i
    enddo
    do i=1,NPHY
       oSI(i)=ow_p(NPHY)+i
    enddo

    do i=1,NPHY
       oLno3(i)=oSI(NPHY)+i
    enddo
    do i=1,NPHY
       oTheta(i)=oLno3(NPHY)+i
    enddo

    do i=1,NPHY
       oQN(i)=oTheta(NPHY)+i
    enddo

    oPPt  =oQN(NPHY)+1
    oD_NO3=oPPt +1

    do i=1,NPHY
       oD_PHY(i)=oD_NO3+i
    enddo

    oD_ZOO=oD_PHY(NPHY)+1
    oD_DET=oD_ZOO+1
    Nout  =oD_DET
    write(6,'(A8,1x, I3, A25)') 'Totally',Nout,'Variables for diagosis!'
    allocate(Varout(  Nout,nlev))
    allocate(Labelout(Nout+ ow ))

    Labelout(oTemp  )='Temp '
    Labelout(oPAR   )='PAR  '
    Labelout(oAks   )='Aks  '
    Labelout(ow     )='w    '
    Labelout(oNO3+ow)='NO3  '

    do i=1,NPHY
       write(Labelout(oPHY(i)+ow), "(A3,I2)") 'PHY',i
    enddo

    Labelout(oZOO +ow)='ZOO  '
    Labelout(oDET +ow)='DET  '
    Labelout(oCHLt+ow)='CHL_T'
    do i=1,NPHY
       write(Labelout(omuNet(i) + ow), "(A3,I2)") 'muN',i
       write(Labelout(oGraz(i)  + ow), "(A3,I2)") 'Gra',i
       write(Labelout(ow_p(i)   + ow), "(A3,I2)") 'w_p',i
       write(Labelout(oD_PHY(i) + ow), "(A3,I2)") 'D_P',i
       write(Labelout(oTheta(i) + ow), "(A3,I2)") 'THE',i
       write(Labelout(oQN(i)    + ow), "(A3,I2)") 'QN ',i
       write(Labelout(oSI(i)    + ow), "(A3,I2)") 'SI ',i
       write(Labelout(oLno3(i)  + ow), "(A3,I2)") 'LNO',i
    enddo
  
    Labelout(oZ2N  + ow)='Z2N  '
    Labelout(oD2N  + ow)='D2N  '
    Labelout(oPPt  + ow)='NPP_T'
    Labelout(oD_NO3+ ow)='D_NO3'
    Labelout(oD_ZOO+ ow)='D_ZOO'
    Labelout(oD_DET+ ow)='D_DET'

    do i = 1, Nout+ow
       write(6,'(A9,I2,A4,A5)') 'Labelout(',i,') = ',Labelout(i)
    enddo
    ! Initialize parameters
    ! Indices for parameters that will be used in MCMC                 
    !imu0    =  1
    !iaI0    =  imu0    + 1
    !iV0N    =  iaI0    + 1
    iaI0 = 1

    if (nutrient_uptake .eq. 1) then
       iKN  =  iaI0    + 1
       iQ0N =  iKN     + 1
    else if(nutrient_uptake .eq. 2) then 
       iA0N =  iaI0    + 1
       iQ0N =  iA0N    + 1
    else
       write(6,*) 'Option of nutrient uptake incorrect! Quit!'
       stop
    endif

    igmax   =  iQ0N    + 1
    ikp     =  igmax   + 1 
    iwDET   =  ikp     + 1
    irdN    =  iwDET   + 1
    imz     =  irdN    + 1
    NPar    =  imz
    !iPHYini =  imz     + 1
    !iZOOini =  iPHYini + 1
    !iDETini =  iZOOini + 1
    !NPar    =  iDETini

    allocate(params(NPar))
    write(6,'(A23,I3,A30)') 'There are totally ',NPar,'parameters to be fitted.'

    allocate(ParamLabel(NPar))
    !ParamLabel(imu0 ) = 'mu0hat '
    ParamLabel(iaI0 ) = 'aI0    '

    if (nutrient_uptake .eq. 1) then
        ParamLabel(iKN  ) = 'KN     '
    else if(nutrient_uptake .eq. 2) then
        ParamLabel(iA0N ) = 'A0N    '
    else
        write(6,*) 'Option of nutrient uptake incorrect! Quit!'
        stop
    endif

    !ParamLabel(iV0N )  = 'V0N    '
    ParamLabel(iQ0N )  = 'Q0N    '
    ParamLabel(iwDET)  = 'wDET   '
    ParamLabel(igmax)  = 'gmax   '
    ParamLabel(ikp  )  = 'kp     '
    ParamLabel(irdn )  = 'rdn    '
    ParamLabel(imz  )  = 'mz     '
    !ParamLabel(iPHYini)= 'PHYini '
    !ParamLabel(iZOOini)= 'ZOOini '
    !ParamLabel(iDETini)= 'DETini '

    !params(imu0)       = 5d0
    params(iaI0)       = 1d0
    !params(iV0N)       = 5d0
    params(igmax)      = 1d0
    params(ikp)        = 1d0
    if (nutrient_uptake .eq. 1) then
       params(iKN)     = 0.1
    elseif (nutrient_uptake .eq. 2) then
       params(iA0N )   = 1d0
    else
       write(6,*) 'Option of nutrient uptake incorrect! Quit!'
       stop
    endif
    params(iQ0N )   = 0.05
    params(iwDET)   = 1d0
    params(irdn )   = 0.05
    params(imz  )   = 0.05
    !params(iPHYini) = 0.1
    !params(iZOOini) = 0.1
    !params(iDETini) = 0.1

  case(EFTdiscrete)

    WRITE(6,*) 'EFTdiscrete model (20 PHY size classes) selected!'
    NPHY=20
    ! Assign variable array indices 
    ! and not including CHL in the variable lists 
    allocate(iPHY(NPHY))

    do i=1,NPHY
       iPHY(i)=i+iNO3
    enddo 

    iZOO=iPHY(NPHY)+1
    iDET=iZOO+1
    NVAR=iDET 

    allocate(Vars(NVAR,nlev),  STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

    NVsinkterms = 1+NPHY
    allocate(Windex(NPHY+1),  STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

    do i=1,NPHY
       Windex(i)=iPHY(i)
    enddo
    Windex(NVsinkterms)=iDET

    ! Output array matrices
    allocate(oPHY(NPHY))
    do i=1,NPHY
       oPHY(i)=i+oNO3
    enddo
    oZOO=oPHY(NPHY)+1
    oDET=oZOO+1
    oFER=oDET+1

    allocate(oTheta(NPHY))
    allocate(oQN(NPHY))
    allocate(omuNet(NPHY))
    allocate(oGraz(NPHY))
    allocate(ow_p(NPHY))
    allocate(oSI(NPHY))
    allocate(oLno3(NPHY))
    allocate(oD_PHY(NPHY))
    allocate(oTheHat(NPHY))

    do i=1,NPHY
       oTheta(i)=oFER+i
    enddo

    do i=1,NPHY
       oQN(i)=oTheta(NPHY)+i
    enddo

    do i=1,NPHY
       omuNet(i)=oQN(NPHY)+i
    enddo

    do i=1,NPHY
       oGraz(i)=omuNet(NPHY)+i
    enddo

    oZ2N=oGraz(NPHY)+1

    oD2N=oZ2N+1
    do i=1,NPHY
       ow_p(i)=oD2N+i
    enddo

    do i=1,NPHY
       oSI(i)=ow_p(NPHY)+i
    enddo

    do i=1,NPHY
       oLno3(i)=oSI(NPHY)+i
    enddo

    do i=1,NPHY
       oTheHat(i)=oLno3(NPHY)+i
    enddo

    oPHYt =oTheHat(NPHY)+1
    oCHLt =oPHYt+1

    do i=1,4
       oCHLs(i)=oCHLt+i
    enddo
    oPPt  =oCHLs(4)+1
    oPMU  =oPPt+1
    oVAR  =oPMU+1
    oD_NO3=oVAR+1

    do i=1,NPHY
       oD_PHY(i)=oD_NO3+i
    enddo

    oD_ZOO=oD_PHY(NPHY)+1
    oD_DET=oD_ZOO+1
    Nout  =oD_DET

    allocate(Varout(    Nout,nlev),  STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

    allocate(Labelout(  Nout+ ow ),  STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

    Labelout(oTemp  )='Temp '
    Labelout(oPAR   )='PAR  '
    Labelout(oAks   )='Aks  '
    Labelout(ow     )='w    '
    Labelout(oNO3+ow)='NO3  '
    do i=1,NPHY
       write(Labelout(oPHY(i)+ow), "(A3,I2)") 'PHY',i
    enddo
    Labelout(oZOO+ow)='ZOO  '
    Labelout(oDET+ow)='DET  '
    Labelout(oFER+ow)='FER  '
    do i=1,NPHY
       write(Labelout(oTheta(i)  + ow), "(A3,I2)") 'THE',i
       write(Labelout(oQN(i)     + ow), "(A3,I2)") 'QN ',i
       write(Labelout(omuNet(i)  + ow), "(A3,I2)") 'muN',i
       write(Labelout(oGraz(i)   + ow), "(A3,I2)") 'Gra',i
       write(Labelout(ow_p(i)    + ow), "(A3,I2)") 'w_p',i
       write(Labelout(oSI(i)     + ow), "(A3,I2)") 'SI ',i
       write(Labelout(oLno3(i)   + ow), "(A3,I2)") 'LNO',i
       write(Labelout(oTheHat(i) + ow), "(A3,I2)") 'THA',i
       write(Labelout(oD_PHY(i)  + ow), "(A3,I2)") 'D_P',i
    enddo
    do i=1,4
       write(Labelout(oCHLs(i)   + ow), "(A4,I2)") 'CHLs',i
    enddo

    Labelout(oZ2N  + ow)='Z2N  '
    Labelout(oD2N  + ow)='D2N  '
    Labelout(oPHYt + ow)='PHY_T'
    Labelout(oCHLt + ow)='CHL_T'
    Labelout(oPPt  + ow)='NPP_T'
    Labelout(oPMU  + ow)='PMU  '
    Labelout(oVAR  + ow)='VAR  '
    Labelout(oD_NO3+ ow)='D_NO3'
    Labelout(oD_ZOO+ ow)='D_ZOO'
    Labelout(oD_DET+ ow)='D_DET'
 
    ! Initialize parameters
    ! Indices for parameters that will be used in MCMC                 
    !imu0    =  1
    !iaI0    =  imu0    + 1
    iaI0    =  1
    ialphaI =  iaI0    + 1
    !iV0N    =  ialphaI + 1
    !ialphaV =  iV0N    + 1
    select case(nutrient_uptake)
    case(1)
      iKN     =  ialphaI + 1
      ialphaK =  iKN     + 1
      iQ0N    =  ialphaK + 1
    case(2)
      iA0N    =  ialphaI + 1
      ialphaA =  iA0N    + 1
      iQ0N    =  ialphaA + 1
    case default
      write(6,*) 'Option of nutrient uptake incorrect! Quit!'
      stop
    end select

    ialphaQ =  iQ0N    + 1
    igmax   =  ialphaQ + 1
    ikp     =  igmax   + 1 
    ialphaG =  ikp     + 1
    iwDET   =  ialphaG + 1
    irdN    =  iwDET   + 1
    imz     =  irdN    + 1
    NPar    =  imz
    !iPHYini =  imz     + 1
    !iZOOini =  iPHYini + 1
    !iDETini =  iZOOini + 1
    !NPar    =  iDETini

    allocate(params(NPar))
    allocate(ParamLabel(NPar))

!    ParamLabel(imu0    ) = 'mu0hat '
    ParamLabel(iaI0    ) = 'aI0    '
    ParamLabel(ialphaI ) = 'alphaI '
    if (nutrient_uptake .eq. 1) then
       ParamLabel(iKN    )  = 'K0N    '
       ParamLabel(ialphaK ) = 'alphaK '
    elseif (nutrient_uptake .eq. 2) then
       ParamLabel(iA0N    ) = 'A0N    '
       ParamLabel(ialphaA ) = 'alphaA '
    else
       print *, 'Option of nutrient uptake incorrect! Quit!'
       stop
    endif
!    ParamLabel(iV0N    ) = 'V0N    '
!    ParamLabel(ialphaV ) = 'alphaV '
    ParamLabel(iQ0N    ) = 'Q0N    '
    ParamLabel(ialphaQ ) = 'alphaQ '
    ParamLabel(ialphaG ) = 'alphaG '
    ParamLabel(iwDET   ) = 'wDET   '
    ParamLabel(igmax   ) = 'gmax   '
    ParamLabel(ikp     ) = 'kp     '
    ParamLabel(irdn    ) = 'rdn    '
    ParamLabel(imz     ) = 'mz     '
    !ParamLabel(iPHYini ) = 'PHYini '
    !ParamLabel(iZOOini ) = 'ZOOini '
    !ParamLabel(iDETini ) = 'DETini '

    !params(iPHYini) =0.1
    !params(iZOOini) =0.1
    !params(iDETini) =0.1
  !  params(imu0)    =5d0
    params(iaI0)    =0.5
  !  params(iV0N)    =5d0
    params(igmax)   =1d0
    params(ikp     )=5d-1
    params(ialphaI) =-0.13
    if (nutrient_uptake .eq. 1) then
       params(iKN)     =1d0
       params(ialphaK) =0.27
    elseif (nutrient_uptake .eq. 2) then
       params(iA0N   ) =1D0
       params(ialphaA) =-0.3d0
    endif
  !  params(ialphaV) =1D-6      ! To avoid zero
    params(ialphaG) =1.1d0
    params(iQ0N   ) =4d-2
    params(ialphaQ) =-0.17d0
    params(iwDET  ) =1D0
    params(irdn   ) =5d-2
    params(imz    ) =5d-2

    call assign_PMU

    CASE(EFTcont)

    NPHY=1
    ! Assign variable array indices 
    ! and not including CHL in the state variable lists 
    allocate(iPHY(NPHY))

    do i=1,NPHY
       iPHY(i)=i+iNO3
    enddo 

    iZOO = iPHY(NPHY)+ 1
    iDET = iZOO      + 1
    iPMU = iDET      + 1
    iVAR = iPMU      + 1
    NVAR = iVAR 

    allocate(Vars(NVAR,nlev))
    NVsinkterms = NPHY+1

    allocate(Windex(NVsinkterms))
    do i=1,NPHY
       Windex(i)=iPHY(i)
    enddo

    Windex(NVsinkterms)=iDET

    ! Output array matrices
    allocate(oPHY(NPHY))
    do i=1,NPHY
       oPHY(i)=i+oNO3
    enddo

    oZOO=oPHY(NPHY)+1
    oDET=oZOO+1
    oPMU=oDET+1
    oVAR=oPMU+1
    oFER=oVAR+1

    allocate(oTheta(NPHY))
    allocate(oQN(NPHY))
    allocate(omuNet(NPHY))
    allocate(oGraz(NPHY))
    allocate(ow_p(NPHY))
    allocate(oD_PHY(NPHY))

    do i=1,NPHY
       oTheta(i)=oFER+i
    enddo

    do i=1,NPHY
       oQN(i)=oTheta(NPHY)+i
    enddo

    do i=1,NPHY
       omuNet(i)=oQN(NPHY)+i
    enddo

    do i=1,NPHY
       oGraz(i)=omuNet(NPHY)+i
    enddo

    oCHLt=oGraz(NPHY)+1
    do i=1,4
       oCHLs(i)=oCHLt+i
    enddo
    oPPt=oCHLs(4)+1
    oZ2N=oPPt+1
    oD2N=oZ2N+1
    do i=1,NPHY
       ow_p(i)=oD2N+i
    enddo

    oD_NO3=ow_p(NPHY)+1
    do i=1,NPHY
       oD_PHY(i)=oD_NO3+1
    enddo
    oD_ZOO=oD_PHY(NPHY)+1
    oD_DET=oD_ZOO+1
    oD_PMU=oD_DET+1
    oD_VAR=oD_PMU+1
    odmudl=oD_VAR+1
    od2mu =odmudl+1
    od2gdl=od2mu +1
    Nout  =od2gdl

    allocate(Varout(Nout,nlev), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate(Labelout(Nout+ow))

    Labelout(oTemp ) ='Temp '
    Labelout(oPAR  ) ='PAR  '
    Labelout(oAks  ) ='Aks  '
    Labelout(ow    ) ='w    '
    Labelout(oNO3+ow)='NO3  '
    do i=1,NPHY
       write(Labelout(oPHY(i)+ow), "(A3,I2)") 'PHY',i
    enddo
    Labelout(oZOO+ow)='ZOO  '
    Labelout(oDET+ow)='DET  '
    Labelout(oPMU+ow)='PMU  '
    Labelout(oVAR+ow)='VAR  '
    Labelout(oFER+ow)='FER  '
    do i=1,NPHY
       write(Labelout(oTheta(i)+ow), "(A3,I2)") 'The',i
       write(Labelout(oQN(i)   +ow), "(A3,I2)") 'QN ',i
       write(Labelout(omuNet(i)+ow), "(A3,I2)") 'muN',i
       write(Labelout(oGraz(i) +ow), "(A3,I2)") 'Gra',i
       write(Labelout(ow_p(i)  +ow), "(A3,I2)") 'w_p',i
       write(Labelout(oD_PHY(i)+ow), "(A3,I2)") 'D_P',i
    enddo
    do i=1,4
       write(Labelout(oCHLs(i) +ow), "(A4,I2)") 'CHLs',i
    enddo

    Labelout(oCHLt +ow)='CHL_T'
    Labelout(oPPt  +ow)='NPP_T'
    Labelout(oZ2N  +ow)='Z2N  '
    Labelout(oD2N  +ow)='D2N  '
    Labelout(oD_NO3+ow)='D_NO3'
    Labelout(oD_ZOO+ow)='D_ZOO'
    Labelout(oD_DET+ow)='D_DET'
    Labelout(oD_PMU+ow)='D_PMU'
    Labelout(oD_VAR+ow)='D_VAR'
    Labelout(odmudl+ow)='dmudl'
    Labelout(od2mu +ow)='d2mu '
    Labelout(od2gdl+ow)='d2gdl'

    do i = 1, Nout+ow
       write(6,*) 'Labelout(',i,') = ',Labelout(i)
    enddo
    ! Initialize parameters
    !===================================================
    ! Define indices for parameters that will be used in MCMC                 
    !===================================================
    !imu0    =  1
    !iaI0    =  imu0    + 1
    iaI0    =  1
    ialphaI =  iaI0    + 1
    !iV0N    =  ialphaI + 1
    !ialphaV =  iV0N    + 1

    select case(nutrient_uptake)
    case(1)
      iKN     =  ialphaI + 1
      ialphaK =  iKN     + 1
      iQ0N    =  ialphaK + 1
    case(2)
      iA0N    =  ialphaI + 1
      ialphaA =  iA0N    + 1
      iQ0N    =  ialphaA + 1
    case default
      write(6,*) 'Option of nutrient uptake incorrect! Quit!'
      stop
    end select

    ialphaQ =  iQ0N    + 1
    igmax   =  ialphaQ + 1
    ikp     =  igmax   + 1 
    ialphaG =  ikp     + 1
    iwDET   =  ialphaG + 1
    irdN    =  iwDET   + 1
    imz     =  irdN    + 1
    iPenfac =  imz     + 1
    ithetamin= iPenfac + 1
    iQNmin  =  ithetamin+1
    NPar    =  iQNmin
    !iPHYini =  imz     + 1
    !iZOOini =  iPHYini + 1
    !iDETini =  iZOOini + 1
    !NPar    =  iDETini

    allocate(params(NPar))
    allocate(ParamLabel(NPar))
    !===================================================
    ! Give parameter labels
    !===================================================

!    ParamLabel(imu0    ) = 'mu0hat '
    ParamLabel(iaI0    ) = 'aI0    '
    ParamLabel(ialphaI ) = 'alphaI '
    if (nutrient_uptake .eq. 1) then
        ParamLabel(iKN     ) = 'K0N    '
        ParamLabel(ialphaK ) = 'alphaK '
    elseif (nutrient_uptake .eq. 2) then
        ParamLabel(iA0N    ) = 'A0N    '
        ParamLabel(ialphaA ) = 'alphaA '
    endif
!    ParamLabel(iV0N    ) = 'V0N    '
!    ParamLabel(ialphaV ) = 'alphaV '
    ParamLabel(iQ0N    ) = 'Q0N    '
    ParamLabel(ialphaQ ) = 'alphaQ '
    ParamLabel(ialphaG ) = 'alphaG '
    ParamLabel(iwDET   ) = 'wDET   '
    ParamLabel(igmax   ) = 'gmax   '
    ParamLabel(ikp     ) = 'kp     '
    ParamLabel(irdn    ) = 'rdn    '
    ParamLabel(imz     ) = 'mz     '
    ParamLabel(iPenfac ) = 'Penfac '
    ParamLabel(ithetamin)= 'thetmin'
    ParamLabel(iQNmin )  = 'QNmin '
    !ParamLabel(iPHYini ) = 'PHYini '
    !ParamLabel(iZOOini ) = 'ZOOini '
    !ParamLabel(iDETini ) = 'DETini '

!    params(imu0    )=5d0
    params(iaI0    )=0.25d0
!    params(iV0N    )=5d0
    params(igmax   )=1d0
    params(ialphaI )=-0.13
    if (nutrient_uptake .eq. 1) then
       params(iKN     )=1d0
       params(ialphaK )=0.27
    elseif (nutrient_uptake .eq. 2) then
       params(iA0N    )=4D1
       params(ialphaA )=-3d-1
    else
       write(6,*) 'Option of nutrient uptake incorrect!'
       stop
    endif
!    params(ialphaV)=2d-1
    params(ialphaG)=1.1d0
    params(iQ0N   )=1d-1
    params(ialphaQ)=-0.17d0
    params(iwDET  )=1D0
    params(ikp    )=1D0
    params(irdn   )=5d-2
    params(imz    )=5d-2
    params(iPenfac)=1d0
    params(ithetamin)=0.05
    params(iQNmin)=1d-2
    !params(iPHYini)=0.1
    !params(iZOOini)=0.1
    !params(iDETini)=0.1

  CASE DEFAULT

    write(6,*) 'Error: Incorrect option for biological models!'
    stop

ENDSELECT
end subroutine choose_model
!========================================================
! Calculate total nitrogen in the system
subroutine Cal_TN
implicit none
integer :: i,k

Ntot = 0d0

do k = 1,nlev
   do i = 1,iDET
      Ntot = Ntot + Vars(i,k) * Hz(k)
   enddo
enddo

end subroutine Cal_TN
!========================================================
subroutine Phygrowth_size(PMU, NO3, tC, par_, muNet, QN, Theta, wPHY, LSI, Lno3, ThetaHat)
implicit none
!INPUT PARAMETERS:
real, intent(in)    :: PMU, NO3, tC, par_
real, intent(out)   :: muNet, QN, Theta, wPHY,  LSI, Lno3, ThetaHat
!LOCAL VARIABLES of phytoplankton:
real :: tf_p,I_zero,X
real :: larg   !Environmental variables
real :: KFe,Fe
real :: V0hat,Kn,A0N,A0hat,fA,muIhat
real :: VNhat
real :: mu0hat,mu0hatSI,RMchl
real :: fV, SI
real :: ZINT
real :: aI
real :: Qs
real :: larg1,w_p,W_Y
logical, parameter :: do_IRON = .false.
!-----------------------------------------------------------------------
!Warning: All rates should be multiplied by dtdays to get the real rate

tf_p    = TEMPBOL(Ep,tC)

! Fe related:

if (do_IRON) then
   Fe   = max(Fe,Femin)  !Dissolved Fe concentration
   KFe  = ScaleTrait(PMU, K0Fe,alphaFe) !The half saturation constant of Fe at average size
endif
    
Qs      = ScaleTrait(PMU, params(iQ0N), params(ialphaQ))/2d0
mu0hat  = dtdays*tf_p*mu0 * exp(alphamu*PMU + betamu*PMU*PMU)

! X: a scratch variable
X       = 0d0

if (do_IRON) then
  mu0hat= mu0hat * Fe/(Fe + KFe)
  X     = alphaFe*KFe/(Fe + KFe) 
endif
  
  ! Iron limits nitrogen uptake
  V0hat=ScaleTrait(PMU, dtdays*tf_p*V0N, alphaV)
  
  if (do_IRON) V0hat = V0hat*Fe/(Fe + KFe)
  
  ! Initial slope of P-I curve
  aI = ScaleTrait(PMU, dtdays*params(iaI0), params(ialphaI))
  ! Cost of photosynthesis
  RMchl  = tf_p*RMchl0*dtdays
  ! Threshold irradiance and RMchl is set temperature dependent
  I_zero = zetaChl*RMchl/aI  
  
  !Define VNhat: the saturation function of ambient nutrient concentration
  SELECTCASE(nutrient_uptake)  
  ! case 1: Classic Michaelis Menton 
    case(1)
  ! Potential maximal nutrient-limited uptake
      ! Half-saturation constant of nitrate uptake
      Kn     = ScaleTrait(PMU,params(iKN), params(ialphaK)) 
      VNhat  = V0hat*NO3/(NO3 + Kn)
  
     ! case 2: optimal uptake based on Pahlow (2005) and Smith et al. (2009)
    case(2)
  
    !  Lmin  = log((params(iLref)**3)/6d0*pi) &
    !  + log(1D0 +params(iPenfac))/( params(iPenfac)*params(ialphaA)) 

      A0N   = dtdays*tf_p*params(iA0N)
      A0hat = ScaleTrait(PMU, A0N, params(ialphaA))

    !  A0hat = PenAff(PMU, params(ialphaA), params(iPenfac),Lmin) &
    !   * ScaleTrait(PMU, A0N, params(ialphaA))

      A0hat = max(A0hat,1D-9*A0N)
 
      !Define fA
      fA = 1D0/( 1D0 + sqrt(A0hat*NO3/V0hat) ) 
    
      VNhat = (1D0-fA)*V0hat*fA*A0hat*NO3/((1D0-fA)*V0hat + fA*A0hat*NO3) 
      
    case default
     write(6,*) 'Error: Incorrect option for nutrient uptake!'
     stop
  ENDSELECT  

! Calculate thetahat (optimal g Chl/mol C for the chloroplast under nutrient replete conditions)
! Only calculate within the euphotic zone, otherwise many numerical problems.

  if(par_ .gt. I_zero) then
    
    larg1 = exp(1d0 + aI*par_/(mu0hat*zetaChl))

    larg  = (1d0 + RMchl/mu0hat)*larg1   
    
   W_Y      = WAPR(larg,0,0)
   ThetaHat = 1d0/zetaChl + (1d0- W_Y)*mu0hat/(aI * par_)

    ! Effect of light limitation
    SI = 1d0 - max(exp(-aI*par_*ThetaHat/mu0hat),0d0)

! Light dependent growth rate 
! (needs to take into account the cost of dark and light-dependent chl maintenance)
    mu0hatSI = mu0hat*SI  ! Gross specific carbon uptake (photosynthesis)
   
    muIhat   = mu0hatSI-(mu0hatSI+RMchl)*zetaChl*ThetaHat ! Net specific carbon uptake
    muIhat   = max(muIhat,1D-9*mu0hat)
   
    LSI    = 1D0-SI

    ZINT   = Qs*(muIhat/VNhat + zetaN)
           
    fV = (-1d0 + sqrt(1d0 + 1d0/ZINT))*Qs*muIhat/VNhat
    fV = max(fV,0.01)
    
    else
    ! Under the conditions of no light:
       ThetaHat      = 1d-2  !  a small positive value 
       ZINT          = Qs*zetaN
       fV            = 1d-2
       muIhat        = 0d0
       LSI           = 1D0
    endif
     ! Nutrient limitation index:
     Lno3 =1d0/(1d0 + sqrt(1D0 +1D0/ZINT)) 

     ! Optimal nutrient quota:
     QN = Qs/Lno3
    
     if (par_ .gt. I_zero) then
!Net growth rate (d-1) of phytoplankton at the average size
      ! muNet = muIhat*(1d0-fV-Qs/QN) - zetaN*fV*VNhat
       muNet = muIhat*(1d0-2d0*Lno3)
     else
       muNet = 0d0
     endif
! Chl:C ratio [ g chl / mol C ] of the whole cell
   Theta  = ThetaHat*(1D0-fV-Qs/QN)
   
! Phytoplankton sinking rate at the average size 
! (should not normalize to dtdays)
   w_p    = ScaleTrait(PMU,abs(w_p0),alphaW)  !Positive
! Constrain the sinking rate not too large for large cells (Smayda 1970)
   wPHY   = min(w_p, Wmax)
end subroutine Phygrowth_size
!========================================================
subroutine assign_PMU
implicit none
integer :: i
real    :: dx
integer :: AllocateStatus
! Calculate the average size 
allocate(PMU_(NPHY), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
PMU_(:)    = 0d0

allocate(wtCHL(4,NPHY), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

wtCHL(:,:) = 0d0

PMU_(1)= PMU_min
dx     = (PMU_max-PMU_min)/float(NPHY-1)

do i=2,NPHY
   PMU_(i) = PMU_(i-1)+dx
enddo

! Calculate the weight for each size class 
! (largest size class being the first to be consistent with observational data)
do i=1,NPHY
   if (PMU_(i) .le. PMU_1) then
      wtCHL(4,i) = 1d0  ! < 1 um
   elseif (PMU_(i) .le. PMU_3) then
      wtCHL(3,i) = 1d0  ! 1-3 um
   elseif (PMU_(i) .le. PMU_10) then
      wtCHL(2,i) = 1d0  ! 3-10 um
   else
      wtCHL(1,i) = 1d0  ! >10 um
   endif
enddo
end subroutine assign_PMU
!========================================================
SUBROUTINE Geider_simple
implicit none
integer :: k
real    :: rhochl, theta, CHL
real    :: rmax_T, SI, muNet, tf_z
real    :: pp_ZP, pp_NZ, pp_ND, pp_DZ
real    :: par_


DO k = 1, nlev

   NO3 = Vars(iNO3,   k)
   PHY = Vars(iPHY(1),k)
   ZOO = Vars(iZOO,   k)
   CHL = Vars(iCHL,   k)
   DET = Vars(iDET,   k)

   ! Check whether in the MLD or not
   if (k .lt. N_MLD) then
      par_ = PAR(k)
   else
      par_ = PARavg
   endif
   ! Loss rate of phytoplankton to detritus equals to zero

   ! Define phytoplankton Chl-to-carbon ratio (gChl/molC)
   theta = CHL/(PHY/params(iQ0N))

   ! The maximal growth rate (rmax_T) under temperature tC 
   rmax_T = params(imu0)* TEMPBOL(Ep,Temp(k)) * dtdays

   !The light limitation index (fpar)
   SI = 1d0 - exp(-params(iaI0) * dtdays *par_ * theta/rmax_T)

   ! Growth rate (d-1) is the product function of temperature, light, and nutrient   
   muNet = rmax_T * (  NO3/(NO3 + params(iKn)) )*SI

   ! define rhochl (gChl/molC;the fraction of phytoplankton carbon production that is devoted to Chl synthesis)
   rhochl = thetm*muNet/(params(iaI0)*dtdays*par_*theta)

   ! The total amount of phytoplankton grazed by zooplankton (molN;gmax is the maximal specific ingestion rate!)
   tf_z   = TEMPBOL(Ez,Temp(k))

   gbar   = grazing(grazing_formulation,params(ikp),PHY)
   INGES  = ZOO*params(igmax)*tf_z * dtdays * gbar

   Zmort = ZOO*ZOO*dtdays* params(imz) *tf_z  !Mortality term for ZOO
 
   !Zooplankton excretion rate (-> DOM)
   RES = INGES*(1d0-GGE-unass)

   !ZOOPLANKTON EGESTION (-> POM)
   EGES = INGES*unass

  ! For production/destruction matrix:
     
  pp_ND = dtdays* params(irDN) *DET*tf_z   
  pp_NZ = ZOO*RES        
  pp_DZ = ZOO*EGES+Zmort 
  pp_ZP = ZOO*INGES      
  
  DET1 = (DET + pp_DZ)/(1d0 + dtdays * params(irDN)*tf_z)
 
  Varout(oDET,k) = max(DET1,1D-9)
  
  Varout(oNO3,k) = (NO3+pp_ND+pp_NZ)/(1d0+ PHY*muNet/NO3)

  Varout(oPHY(1),k) = PHY*(1d0 + muNet)/(1d0 + pp_ZP/PHY)

  Varout(omuNet(1), k) = muNet/dtdays

  Varout(oPPt,k) = PHY*muNet/params(iQ0N)/dtdays

  Varout(oGraz(1), k) = pp_ZP/PHY/dtdays

  Varout(oZOO,k) = (ZOO+pp_ZP)/(1d0+EGES+ZOO*dtdays*params(imz)*tf_z+RES)
    
  Varout(oZ2N,k) = pp_NZ/dtdays
  Varout(oD2N,k) = pp_ND/dtdays

  Varout(oCHLt,k) = (CHL + rhochl*muNet*PHY/params(iQ0N)) &
                 / (1d0 + pp_ZP/PHY)

  Varout(ow_p(1),k) = w_p0   !Should not be converted by dtdays
  
enddo
End subroutine Geider_simple
!========================================================
SUBROUTINE FlexEFT_simple
implicit none
integer :: k
real    :: Qs,QN, mu0hat,V0hat, aI, RMchl
real    :: tf_p,tf_z,I_zero,VNhat
real    :: A0hat,fA,larg1,larg
real    :: W_Y,ThetaHat,SI,mu0hatSI,muIhat,par_
real    :: Lno3
real    :: pp_ND,pp_NZ,pp_PN,pp_DZ,pp_ZP
real    :: ZINT,fV,muNet,gbar

DO k=1,nlev
   DET  = Vars(iDET,k)
   NO3  = Vars(iNO3,k)
   PHY  = Vars(iPHY(1),k)
   ZOO  = Vars(iZOO,k)  !Zooplankton biomass

   ! Check whether in the MLD or not
   if (k .lt. N_MLD) then
      par_ = PAR(k)
   else
      par_ = PARavg
   endif

! Phytoplankton section:
   tf_p = TEMPBOL(Ep,Temp(k))
   Qs   = params(iQ0N)/2d0
 mu0hat = dtdays*tf_p*mu0
  V0hat = dtdays*tf_p*V0N
  ! Initial slope of P-I curve
  aI    = dtdays*params(iaI0)
  ! Cost of photosynthesis
  RMchl  = tf_p*RMchl0*dtdays
  ! Threshold irradiance and RMchl is set temperature dependent
  I_zero = zetaChl*RMchl/aI  
  
!Define VNhat: the saturation function of ambient nutrient concentration
  SELECTCASE(nutrient_uptake)  
  ! case 1: Classic Michaelis Menton 
    case(1)
  ! Potential maximal nutrient-limited uptake
      ! Half-saturation constant of nitrate uptake
      VNhat  = V0hat*NO3/(NO3 + params(iKN))

 ! case 2: optimal uptake based on Pahlow (2005) and Smith et al. (2009)
    case(2)
      A0hat = dtdays*tf_p*params(iA0N)
      !Define fA
      fA = 1D0/(1D0 + sqrt(A0hat*NO3/V0hat) ) 
      VNhat = (1D0-fA)*V0hat*fA*A0hat*NO3/((1D0-fA)*V0hat+fA*A0hat*NO3) 
    case default
     write(6,*) 'Error: Incorrect option for nutrient uptake!'
     stop
  ENDSELECT  

! Calculate thetahat (optimal g Chl/mol C for the chloroplast under nutrient replete conditions)
! Only calculate within the euphotic zone, otherwise many numerical problems.
  if( par_ .gt. I_zero ) then
    larg1 = exp(1d0 + aI*par_/(mu0hat*zetaChl))
    larg  = (1d0 + RMchl/mu0hat)*larg1   
    
     W_Y  = WAPR(larg,0,0)

   ThetaHat = 1d0/zetaChl + (1d0- W_Y)*mu0hat/(aI * par_)

    ! Effect of light limitation
    SI = 1d0 - max(exp(-aI * par_ * ThetaHat/mu0hat),0d0)
    Varout(oSI(1),k) = SI
! Light dependent growth rate 
! (needs to take into account the cost of dark and light-dependent chl maintenance)
    mu0hatSI = mu0hat*SI !Gross specific carbon uptake (photosynthesis)
   
    muIhat   = mu0hatSI-(mu0hatSI+RMchl)*zetaChl*ThetaHat ! Net specific carbon uptake
    muIhat   = max(muIhat,0d0)
   
    ZINT   = Qs*(muIhat/VNhat + zetaN)
           
    fV = (-1d0 + sqrt(1d0 + 1d0/ZINT))*Qs*muIhat/VNhat
    fV = max(fV,0.01)
    
    else
    ! Under the conditions of no light:
       ThetaHat      = 1d-2  !  a small positive value 
       ZINT          = Qs*zetaN
       fV            = 1d-2
       muIhat        = 0d0
       Varout(oSI(1),k) = 0d0
    endif

     ! Nutrient limitation index:
     Lno3 =1d0/(1d0 + sqrt(1D0 +1D0/ZINT)) 
     Varout(oLno3(1),k) = Lno3
     ! Optimal nutrient quota:
     QN   = Qs/Lno3
     Varout(oQN(1),k)   = QN
    
     if (par_ .gt. I_zero) then
   !Net growth rate (d-1) of phytoplankton at the average size
       muNet = muIhat*(1d0-2d0*Lno3)  !Droop model
     else
       muNet = 0d0
     endif

   !Save net growth rate
    Varout(omuNet(1),k) = muNet/dtdays

   !Chl:C ratio [ g chl / mol C ] of the whole cell
    Varout(oTheta(1),k) = ThetaHat*(1D0-fV-Qs/QN)
    
   !Phytoplankton sinking rate at the average size
    Varout(ow_p(1),  k) = abs(w_p0)  !Positive

!---------------------------------------------------------------
!! ZOOplankton section:
   tf_z = TEMPBOL(Ez,Temp(k))

 ! The grazing dependence on total prey (dimensionless)
   gbar = grazing(grazing_formulation,params(ikp),PHY)
 !Zooplankton total ingestion rate
   INGES = tf_z*dtdays*params(igmax)*gbar

 !Zooplankton excretion rate (-> DOM)
   RES  = INGES*(1d0-GGE-unass)

 !ZOOPLANKTON EGESTION (-> POM)
   EGES = INGES*unass
    
! Grazing rate on PHY(specific to N-based Phy biomass, unit: d-1)

   ! Calculate the specific grazing rate for PHY
   Varout(oGraz(1),k) = INGES*ZOO/PHY/dtdays
   Varout(oPHY(1),k)  = PHY*(1d0+muNet)/(1d0 + INGES*ZOO/PHY)
!!End of zooplankton section
!=============================================================
!! Solve ODE functions:
  Zmort = ZOO*ZOO*dtdays*params(imz)*tf_z  !Mortality term for ZOO
 
  ! For production/destruction matrix:
  pp_PN=PHY*muNet
  pp_ND=dtdays* params(irDN)*DET*tf_z   
  pp_NZ=ZOO*RES        
  pp_DZ=ZOO*EGES+Zmort 
  pp_ZP=ZOO*INGES      
  

  Varout(oDET,k) = (DET+pp_DZ)/(1d0 + dtdays * params(irDN)*tf_z)
  Varout(oNO3,k) = (NO3+pp_ND+pp_NZ)/(1d0+pp_PN/NO3)

  Varout(oZOO,k) = (ZOO+pp_ZP)/(1d0   &
                 + EGES+ZOO*dtdays*params(imz)*tf_z+RES)
    
  Varout(oZ2N,k) = pp_NZ/dtdays
  Varout(oD2N,k) = pp_ND/dtdays
  Varout(oPPt,k) = pp_PN/dtdays/QN
  Varout(oCHLt,k)= Vars(iPHY(1),k)/Varout(oQN(1),k)*Varout(oTheta(1),k)
ENDDO
End subroutine FlexEFT_simple
!------------------------
SUBROUTINE FLEXEFT_DISC
implicit none
integer         :: i,j,k
real :: dx,tf_z,par_
real :: PHYtot = 0d0
real :: PHYtot2= 0d0
real :: pp_ND,pp_NZ, pp_DZ,pp_ZP 
real :: ppC_PN = 0d0
real :: pp_PN  = 0d0
real :: P_PMU  = 0d0
real :: P_VAR  = 0d0
real :: CHLtot = 0d0
real :: CHLs(4)= 0d0  ! Size fractionated Chl

DO k=1,nlev
   PHYtot  =0d0  ! Calculate total PHY biomass
   PHYtot2 =0d0  ! total P**alphaG
   pp_PN   =0d0  ! total primary production (mmol N per d per m3)
   ppC_PN  =0d0  ! total primary production (mmol C per d per m3)
   P_PMU   =0d0  ! biomass* logsize
   P_VAR   =0d0  ! biomass* logsize**2
   CHLtot  =0d0  ! Total CHL A
   CHLs(:) =0d0  ! Size fractionated CHL
  
   DET  = Vars(iDET,k)
   ZOO  = Vars(iZOO,k)  !Zooplankton biomass
   NO3  = Vars(iNO3,k)

   ! Check whether in the MLD or not
   if (k .lt. N_MLD) then
      par_ = PAR(k)
   else
      par_ = PARavg
   endif

!! Phytoplankton section:
   do i=1,NPHY

     call Phygrowth_size(PMU_(i),Vars(iNO3,k),Temp(k),par_, &
          Varout(omuNet(i),k), Varout(oQN(i),k), Varout(oTheta(i),k),&
          Varout(ow_p(i)  ,k), Varout(oSI(i),k), Varout(oLno3(i), k),&
          Varout(oTheHat(i),k))

     PHYtot =PHYtot +Vars(iPHY(i),k)

     PHYtot2=PHYtot2+Vars(iPHY(i),k)**params(ialphaG)

     pp_PN  =Vars(iPHY(i),k)*Varout(omuNet(i),k)+pp_PN

     ! Chl in this size class
     dx     =Vars(iPHY(i),k)/Varout(oQN(i),k)   *Varout(oTheta(i),k)
     
     ppC_PN =Vars(iPHY(i),k)*Varout(omuNet(i),k)/Varout(oQN(i),k)+ppC_PN
     CHLtot =CHLtot+ dx
     P_PMU  =P_PMU + Vars(iPHY(i),k)*PMU_(i)
     P_VAR  =P_VAR + Vars(iPHY(i),k)*PMU_(i)**2

     do j = 1, 4
        CHLs(j) = CHLs(j) + dx*wtCHL(j,i)
     enddo
   enddo

   ! save mean size
   Varout(oPMU,k) = P_PMU/PHYtot        

   ! save size variance
   Varout(oVAR,k) = P_VAR/PHYtot - Varout(oPMU,k)**2

   ! save total phytoplankton biomass
   Varout(oPHYt,k) = PHYtot
   Varout(oCHLt,k) = CHLtot
   ! save total NPP (carbon-based)
   Varout(oPPt,k)  = ppC_PN/dtdays

   ! save size fractionated Chl
   do j = 1,4
      Varout(oCHLs(j),k) = CHLs(j)
   enddo
!---------------------------------------------------------------
!! ZOOplankton section:
   tf_z = TEMPBOL(Ez,Temp(k))

 ! The grazing dependence on total prey (dimensionless)
   gbar = grazing(grazing_formulation,params(ikp),PHYtot)
 !Zooplankton total ingestion rate
   INGES = tf_z*dtdays*params(igmax)*gbar

 !Zooplankton excretion rate (-> DOM)
   RES  = INGES*(1d0-GGE-unass)

 !ZOOPLANKTON EGESTION (-> POM)
   EGES = INGES*unass
    
! Grazing rate on PHY each size class (specific to N-based Phy biomass, unit: d-1) (Eq. 12)

   ! Calculate the specific grazing rate for each size class
   do i=1,NPHY
    
     ! Eq. 10 in Smith & Sergio
     Varout(oGraz(i),k) = (INGES*ZOO/PHYtot2)   &
       *Vars(iPHY(i),k)**(params(ialphaG)-1d0)
     
     Varout(oPHY(i),k) = Vars(iPHY(i),k)*(1d0+Varout(omuNet(i),k))  &
       /(1d0 + Varout(oGraz(i),k))
      
   enddo

!!End of zooplankton section
!=============================================================
!! Solve ODE functions:
  Zmort = ZOO*ZOO*dtdays* params(imz) *tf_z  !Mortality term for ZOO
 
  ! For production/destruction matrix:
  pp_ND = dtdays* params(irDN) *DET*tf_z   
  pp_NZ = ZOO*RES        
  pp_DZ = ZOO*EGES+Zmort 
  pp_ZP = ZOO*INGES      
  DET1  = (DET+pp_DZ)/(1d0 + dtdays * params(irDN)*tf_z)
 
  Varout(oDET,k) = DET1
  Varout(oNO3,k) = (NO3+pp_ND+pp_NZ)/(1d0+pp_PN/NO3)

  Varout(oZOO,k) = (ZOO+pp_ZP)/(1d0   &
                 + EGES+ZOO*dtdays*params(imz)*tf_z+RES)
    
  Varout(oZ2N,k) = pp_NZ/dtdays
  Varout(oD2N,k) = pp_ND/dtdays
ENDDO
END SUBROUTINE FLEXEFT_DISC
!========================================================
SUBROUTINE FLEXEFT_CONT
implicit none
!INPUT PARAMETERS:
real :: tC,par_
!LOCAL VARIABLES of phytoplankton:
real :: PMUPHY,VARPHY,NO31,PHY1,ZOO1
real :: PMUPHY1,VARPHY1  !State variables
real :: larg 
real :: PMU,VAR,PMU1,VAR1,X,B,dXdl,dBdl,KFe,Fe,dmu0hatdl, d2mu0hatdl2
real :: V0hat,Kn,Lmin,A0N,A0N0, A0hat,dA0hatdl,d2A0hatdl2,fA
real :: VNhat,dVNhatdl,d2VNhatdl2 ! Nutrient uptake variables
real :: mu0hat,muIhat,mu0hatSI,dmu0hatSIdl,d2mu0hatSIdl2
real :: dmuIhatdl,d2muIhatdl2! Growth rate variables
real :: fV,dfVdl,d2fVdl2  !fV
real :: ZINT,dZINdl,d2ZINdl2 !ZINT
real :: aI,SI,dSIdl,d2SIdl2         !Light dependent component
real :: RMchl,Theta,ThetaHat,dThetaHatdl,d2ThetaHatdl2 !Chl related variables
real :: QN,Qs,dQNdl,d2QNdl2  ! cell quota related variables
real :: larg1,w_p,dmu0hat_aIdl,d2mu0hat_aIdl2,dlargdl
real :: d2largdl2,W_Y,dWYYdl,daI_mu0hatdl,d2aI_mu0hatdl2
real :: dV0hatdl,d2V0hatdl2,muNet,dmuNetdl,d2muNetdl2
real :: ThetaAvg,QNavg,dwdl, d2wdl2,w_pAvg,dgdlbar,d2gdl2bar, alphaA
real :: Q0N,alphaQ,aI0,alphaI,Lref,Penfac,tf_p,tf_z,I_zero
real :: alphaK, K0N
!Declarations of zooplankton:
integer :: k,j

real :: pp_ND,pp_NZ,pp_PN,pp_NP,pp_DZ,pp_ZP, PPpn
real :: CHLtot  = 0d0
real :: pCHL(4) = 0d0
logical, parameter :: do_IRON = .false.

!-----------------------------------------------------------------------

DO k=1,nlev   
   ! Retrieve current (local) state variable values.
   tC     = Temp(k)
   ! Check whether in the MLD or not
   if (k .lt. N_MLD) then
      par_ = PAR(k)
   else
      par_ = PARavg
   endif

   NO3    = Vars(iNO3,k)
   PHY    = Vars(iPHY(1),k)
   ZOO    = Vars(iZOO,k)
   DET    = Vars(iDET,k)
   PMUPHY = Vars(iPMU,k)
   VARPHY = Vars(iVAR,k)
!All rates have been multiplied by dtdays to get the real rate correponding to the actual time step
   tf_p   = TEMPBOL(Ep,tC)
   DET    = max(DET,1D-6)
   PHY    = max(PHY,1D-6)
   PMUPHY = max(PMUPHY,1D-6)
   VARPHY = max(VARPHY,1D-6)
   PMU    = PMUPHY/PHY - log(1d1)  ! Correct PMU to the real one
   VAR    = VARPHY/PHY
! Fe related:
   if (do_IRON) then
!      Fe   = Vars(iFER,k)

      !Dissolved Fe concentration
      Fe   = max(Fe,Femin)  

      !The half saturation constant of Fe at average size
      KFe  = ScaleTrait(PMU, K0Fe,alphaFe)
   endif
       
   alphaQ  = params(ialphaQ)
   Q0N     = params(iQ0N)
   Qs      = ScaleTrait(PMU, Q0N, alphaQ)/2d0

   mu0hat  = dtdays*tf_p*mu0*exp(alphamu*PMU + betamu*PMU**2)

   X       = 0d0
   
   if (do_IRON) then
     mu0hat  = mu0hat * Fe/(Fe + KFe)
     X       = alphaFe*KFe/(Fe + KFe) 
   endif
   
   dmu0hatdl  = mu0hat*(alphamu + 2D0*betamu*PMU - X)
   d2mu0hatdl2= dmu0hatdl*(alphamu + 2D0*betamu*PMU-X) + mu0hat*2D0*betamu
   
   if (do_IRON) d2mu0hatdl2 = d2mu0hatdl2 - mu0hat*alphaFe*Fe*X/(Fe+KFe)
   
   ! Iron limits nitrogen uptake

   !V0N    = params(iV0N)
   !alphaV = params(ialphaV)
   V0hat  = ScaleTrait(PMU, dtdays*tf_p*V0N, alphaV)
     
   if (do_IRON) V0hat = V0hat*Fe/(Fe + KFe)
     
   dV0hatdl   = V0hat*(alphaV- X)
   d2V0hatdl2 = dV0hatdl*(alphaV - X)
     
   if (do_Iron) d2V0hatdl2 = d2V0hatdl2 - V0hat*alphaFe*Fe*X/(Fe+KFe) 
  
   ! Initial slope of P-I curve
   aI0    = params(iaI0)
   alphaI = params(ialphaI)
   aI     = ScaleTrait(PMU, dtdays*aI0, alphaI)

   ! Cost of photosynthesis
   RMchl  = tf_p*RMchl0*dtdays

   ! Threshold irradiance and RMchl is set temperature dependent
   I_zero = zetaChl*RMchl/aI  
  
   !Define VNhat: the saturation function of ambient nutrient concentration
   SELECTCASE(nutrient_uptake)  
   ! case 1: Classic Michaelis Menton 
   case(1)
   ! Potential maximal nutrient-limited uptake
   ! Half-saturation constant of nitrate uptake
    alphaK = params(ialphaK)
    K0N    = params(iKN)
    Kn     = ScaleTrait(PMU,K0N, alphaK) 
    VNhat  = V0hat*NO3/(NO3 + Kn)
   
    dVNhatdl  = -VNhat*alphaK*Kn/(NO3+Kn) + NO3/(NO3 + Kn)*dV0hatdl
   
     d2VNhatdl2 = -alphaK*(VNhat * alphaK * NO3/(NO3+Kn)**2*Kn          &
   + Kn/(NO3 + Kn)*dVNhatdl)+ NO3/(NO3+ Kn)*d2V0hatdl2                  &
   - dV0hatdl*NO3/(NO3 + Kn)**2*alphaK*Kn
  
     ! case 2: optimal uptake based on Pahlow (2005) and Smith et al. (2009)
   case(2)

     Lref  = 0.5       !params(iLref) 
     Penfac= params(iPenfac)
     alphaA= params(ialphaA)

     Lmin  = log(Lref**3/6d0*pi) + log(1d0+Penfac)/( Penfac*alphaA) 

     A0N0  = params(iA0N)
     A0N   = dtdays*tf_p*A0N0
     A0hat = PenAff(PMU, alphaA, Penfac,Lmin)*ScaleTrait(PMU,A0N,alphaA)
     A0hat = max(A0hat,1D-16*A0N)  ! Maintain positivity
  
     dA0hatdl = alphaA*A0hat-A0N*exp(PMU*alphaA)*Penfac*alphaA           &
             * exp(Penfac*alphaA *(PMU-Lmin))
    
     d2A0hatdl2 = alphaA*dA0hatdl                                        &
             - Penfac*alphaA*exp(alphaA*((1d0+Penfac)*PMU-Penfac*Lmin))  &
             * (dA0hatdl + A0hat*alphaA*(1d0+Penfac))  
      
             !Define fA
      fA    = 1D0/( 1D0 + sqrt(A0hat * NO3/V0hat) ) 
    
      VNhat = (1D0-fA)*V0hat*fA*A0hat*NO3/((1D0-fA)*V0hat+fA*A0hat*NO3) 
      
      !X: temporary variable
      X    = V0hat/A0hat + 2d0*sqrt(V0hat*NO3/A0hat) + NO3       
    
      !B: d(V0/A0)dl
      B    = dV0hatdl/A0hat - V0hat/A0hat**2*dA0hatdl
    
      dXdl = B*(1d0 + sqrt(NO3*A0hat/V0hat))
    
      dBdl = d2V0hatdl2/A0hat - dV0hatdl*dA0hatdl/A0hat**2              &
       - (V0hat/A0hat**2*d2A0hatdl2                                     &
       +  dA0hatdl*(dV0hatdl/A0hat**2 - 2d0*V0hat*dA0hatdl/A0hat**3))
      
      dVNhatdl = NO3*(dV0hatdl/X-V0hat/X**2*B*(1d0+ sqrt(NO3*A0hat/V0hat) ) )
    
      d2VNhatdl2 = NO3*(d2V0hatdl2/X - dV0hatdl*dXdl/X**2               &
       - (V0hat/X**2*B*(-sqrt(NO3)/2d0*(A0hat/V0hat)                    &
       * sqrt(A0hat/V0hat) * B)                                         &
       + B*(1d0 + sqrt(NO3*A0hat/V0hat) )                               &
       * (dV0hatdl/X**2 - 2d0*V0hat*dXdl/X**3)                          &
       + V0hat/X**2*(1d0+sqrt(NO3*A0hat/V0hat))*dBdl))
    
   CASE DEFAULT
      write(6,*) 'Error: Incorrect option for nutrient uptake!'
      STOP
   ENDSELECT  

! Calculate thetahat (optimal g Chl/mol C for the chloroplast under nutrient replete conditions)
! Only calculate within the euphotic zone, otherwise many numerical problems.

      if( par_ .gt. I_zero ) then
        
        larg1 = exp(1d0 + min(aI*par_/(mu0hat*zetaChl),6d2))
   
        larg  = (1d0 + RMchl/mu0hat)*larg1   
        
        dmu0hat_aIdl   = (dmu0hatdl - alphaI*mu0hat)/aI
   
        d2mu0hat_aIdl2 = d2mu0hatdl2/aI -alphaI/aI*dmu0hatdl-aI*dmu0hat_aIdl
   
        daI_mu0hatdl = -(aI/mu0hat)**2*dmu0hat_aIdl
   
        d2aI_mu0hatdl2 = -((aI/mu0hat)**2*d2mu0hat_aIdl2                 &
       - 2d0/(mu0hat/aI)**3*(dmu0hat_aIdl)**2)
   
        dlargdl = -RMchl*larg1/(mu0hat**2)*dmu0hatdl                     &
       + (1d0+RMchl/mu0hat)*larg1 * par_/zetaChl*daI_mu0hatdl
        
        d2largdl2 = -RMchl*(larg1*mu0hat**(-2)*d2mu0hatdl2               &
       + larg1*par_/zetaChl*daI_mu0hatdl*mu0hat**(-2)*dmu0hatdl          &
       + larg1*dmu0hatdl*(-2d0*mu0hat**(-3)*dmu0hatdl))                  &
       + par_/zetaChl*((1.+RMchl/mu0hat)*larg1*d2aI_mu0hatdl2            &
       + (1.+RMchl/mu0hat)*larg1*par_/zetaChl*daI_mu0hatdl*daI_mu0hatdl  &
       + RMchl*(-mu0hat**(-2)*dmu0hatdl)*larg1*daI_mu0hatdl)
   
       W_Y      = WAPR(larg,0,0)
       ThetaHat = 1d0/zetaChl + (1d0- W_Y)*mu0hat/(aI * par_)
       ThetaHat = max(ThetaHat,0.01) 
   
       dThetaHatdl = 1d0/par_*(-W_Y/larg/(1.+W_Y)*dlargdl*mu0hat/aI      &
     +  (1d0-W_Y)*dmu0hat_aIdl)
       
       dWYYdl = dlargdl*(-W_Y**2/larg**2/(1d0+W_Y)**3                    &
       -  W_Y/larg**2/(1d0+W_Y) + W_Y/(larg*(1d0+W_Y))**2)
   
       d2ThetaHatdl2 = 1d0/par_*(-(W_Y/larg/(1d0+W_Y)*dlargdl*dmu0hat_aIdl  &
       +  W_Y/larg/(1d0+W_Y)*d2largdl2*mu0hat/aI   &
       +  dWYYdl*dlargdl*mu0hat/aI)               &
       -  W_Y/larg/(1d0+W_Y)*dlargdl * dmu0hat_aIdl &
       +  (1d0-W_Y)*d2mu0hat_aIdl2)
   
        SI = 1d0 - max(exp(-aI*par_*ThetaHat/mu0hat),0d0)
   
        dSIdl = ( (alphaI- dmu0hatdl/mu0hat)   &
        * ThetaHat + dThetaHatdl) * (1d0-SI)*aI*par_/mu0hat    !confirmed
   
       d2SIdl2 = par_*(- dSIdl*aI/mu0hat*(ThetaHat*alphaI     &
      - ThetaHat/mu0hat*dmu0hatdl + dThetaHatdl) + (1d0-SI)   &
      * (ThetaHat*alphaI- ThetaHat/mu0hat*dmu0hatdl + dThetaHatdl) &
      * daI_mu0hatdl + (1d0-SI)*aI/mu0hat*(                 &
      - (d2mu0hatdl2/mu0hat - dmu0hatdl**2/mu0hat**2)*ThetaHat  &
      + (alphaI-dmu0hatdl/mu0hat)*dThetaHatdl + d2ThetaHatdl2)  )
   
       ! Light dependent growth rate 
       ! (needs to take into account the cost of dark and light-dependent chl maintenance)
       mu0hatSI = mu0hat*SI  ! Gross specific carbon uptake (photosynthesis)
   
       muIhat   = mu0hatSI-(mu0hatSI+RMchl)*zetaChl*ThetaHat ! Net specific carbon uptake
       muIhat   = max(muIhat,1D-16*mu0hatSI)
   
       dmu0hatSIdl = SI*dmu0hatdl + mu0hat*dSIdl
   
       d2mu0hatSIdl2 = d2mu0hatdl2*SI+2d0*dmu0hatdl*dSIdl+mu0hat*d2SIdl2 !Correct
   
       dmuIhatdl = (1D0-zetaChl*ThetaHat)*dmu0hatSIdl - dThetaHatdl*zetaChl*(mu0hatSI+RMchl) !Correct
    
       d2muIhatdl2=d2mu0hatSIdl2 &
       -zetaChl*(ThetaHat*d2mu0hatSIdl2+2d0*dThetaHatdl*dmu0hatSIdl   &
       +mu0hatSI*d2ThetaHatdl2)-zetaChl*RMchl*d2ThetaHatdl2   !Correct
    
       ZINT   = Qs*(muIhat/VNhat + zetaN)
    
       dZINdl = Qs*(dmuIhatdl/VNhat - muIhat*dVNhatdl/VNhat**2)+alphaQ*ZINT    
       
       d2ZINdl2 = Qs/VNhat*((alphaQ-dVNhatdl/VNhat)*dmuIhatdl & 
       + d2muIhatdl2) - Qs/VNhat**2*(muIhat*d2VNhatdl2        &
       + dVNhatdl*(dmuIhatdl+alphaQ*muIhat-2d0*muIhat          &
       / VNhat*dVNhatdl)) + alphaQ*dZINdl
       
       fV = (-1d0 + sqrt(1d0 + 1d0/ZINT))*Qs*muIhat/VNhat
       fV = max(fV, 0.01)
    !
       else
    ! Under the conditions of no light:
          ThetaHat      = 0.01  !  a small positive value 
          dThetaHatdl   = 0d0
          d2ThetaHatdl2 = 0d0
          ZINT          = Qs*zetaN
          dZINdl        = alphaQ*ZINT
          d2ZINdl2      = alphaQ*dZINdl
          fV            = 0.01
          muIhat        = 0d0
          dmuIhatdl     = 0d0
          d2muIhatdl2   = 0d0
       endif
    
          ! Optimal nutrient quota:
       QN = (1d0+ sqrt(1d0+1d0/ZINT))*Qs
    
       dQNdl  = alphaQ*QN-dZINdl*Qs/(2d0*ZINT*sqrt(ZINT*(1d0+ZINT))) !confirmed  
    
       d2QNdl2 = alphaQ*dQNdl - Qs/(2d0*ZINT*sqrt(ZINT*(ZINT+1d0)))  &
        *(d2ZINdl2+alphaQ*dZINdl-(2d0*ZINT+1.5d0)/(ZINT*(ZINT+1d0))*dZINdl**2)      ! Confirmed
    
       dfVdl = alphaQ*Qs*(1d0/QN+2d0*zetaN)-(zetaN+Qs/QN**2)*dQNdl  !Confirmed
    
       d2fVdl2 = (alphaQ**2)*Qs*(1d0/QN + 2d0*zetaN)                    &
         -  2d0*alphaQ*Qs*dQNdl/QN**2 + 2d0*(dQNdl**2)*Qs/QN**3         &     
         -     (zetaN + Qs/QN**2) * d2QNdl2  ! Confirmed
    
       if (par_ .gt. I_zero) then

          X=1d0-fV-Qs/QN

        !Net growth rate (d-1) of phytoplankton at the average size
          muNet = muIhat*X - zetaN*fV*VNhat

!Here the derivative of muNet includes respiratory costs of both N Assim and Chl maintenance       
          dmuNetdl = muIhat*(Qs/(QN**2)*dQNdl                           &
      -  alphaQ*Qs/QN-dfVdl) +  X*dmuIhatdl                             & 
      -  zetaN*(fV*dVNhatdl+VNhat*dfVdl)

         d2muNetdl2 = (Qs/(QN*QN))*dQNdl*dmuIhatdl                      &
     + muIhat*(Qs/QN**2)*(dQNdl*(alphaQ - 2d0*dQNdl/QN) + d2QNdl2)      & 
     - (alphaQ*(Qs/QN*(dmuIhatdl + alphaQ*muIhat)                       &
     - muIhat*Qs/(QN**2)*dQNdl))                                        &
     - (muIhat*d2fVdl2+dfVdl*dmuIhatdl)                                 &
     + dmuIhatdl*(Qs/(QN**2)*dQNdl-alphaQ*Qs/QN-dfVdl)  & 
     + X*d2muIhatdl2                   &
     - zetaN*(fV*d2VNhatdl2 + 2d0*dfVdl*dVNhatdl + VNhat*d2fVdl2)  !dC

       else
        muNet      = 0d0
        dmuNetdl   = 0d0
        d2muNetdl2 = 0d0
       endif
!  chl:C ratio [ g chl / mol C ] of the whole cell, at the mean cell size 
       Theta = ThetaHat*X
! Calculate the mean chl:C ratio of phytoplankton (averaged over all sizes) 
    ! using the second derivative of the cell quota w.r.t. log size. (Why negative?)
       ThetaAvg = max(Theta + (VAR/2d0)*(d2ThetaHatdl2*X                &     
       - ThetaHat*(d2fVdl2 + Qs*(dQNdl**2*(2d0/QN)-d2QNdl2)/(QN**2))),  &
       params(ithetamin))
                
    ! Calculate the mean N:C ratio of phytoplankton (averaged over all sizes) 
    ! using the second derivative of the cell quota w.r.t. log size. (How to derive?)
     QNAvg = max(QN/(1d0+((2d0/QN)*dQNdl**2 - d2QNdl2)*VAR/(2d0*QN)),     &
        params(iQNmin))  
!=============================================================
    ! Calculate the community sinking rate of phytoplankton
    ! Phytoplankton sinking rate at the average size
    w_p    = ScaleTrait(PMU,abs(w_p0),alphaW)  !Positive

    ! Constrain the sinking rate not too large for large cells (Smayda 1970)
    w_p    = min(w_p, Wmax)
    dwdl   = alphaW*w_p
    d2wdl2 = alphaW*dwdl
 
  ! Sinking rate (m) of phytoplankton
    w_pAvg = w_p+0.5*VAR*d2wdl2
  ! sinking rate must be negative!
    w_pAvg = -min(w_pAvg, Wmax)
    
!=============================================================
!! ZOOplankton section:
   tf_z = TEMPBOL(Ez,tC)
  
!   ! The grazing dependence on total prey
   gbar = grazing(grazing_formulation,params(ikp),PHY)

   !Zooplankton per capita total ingestion rate
    INGES = tf_z*dtdays*params(igmax)*gbar

   !Zooplankton excretion rate (-> DOM)
   RES = INGES*(1d0-GGE-unass)

   !ZOOPLANKTON EGESTION (-> POM)
   EGES = INGES*unass
    
! Grazing rate on the mean size of PHY (specific to N-based Phy biomass, unit: d-1) (Eq. 12)
   gbar      = INGES*ZOO/PHY*sqrt(params(ialphaG))
   dgdlbar   = 0d0
   d2gdl2bar = (1d0-params(ialphaG))/VAR*gbar
!!End of zooplankton section
!=============================================================
!! Solve ODE functions:
   Zmort = ZOO*ZOO*dtdays*params(imz)*tf_z  !Mortality term for ZOO

   if (d2gdl2bar .lt. 0d0) then
      VAR1 = VAR*(1d0-VAR*d2gdl2bar)
   else                                       
      VAR1 = VAR/(1d0+VAR*d2gdl2bar)
   endif

   ! Update PMU and VAR:    
   if (d2muNetdl2 .gt. 0d0) then
! self%dVARdt=VAR*VAR*(self%d2muNetdl2- self%d2gdl2bar-self%d2wdl2)  !Eq. 19 
      VAR1 =VAR1*(1d0+VAR*d2muNetdl2)
   else    
      VAR1 =VAR1/(1d0-VAR*d2muNetdl2)
   endif
   VAR1    =VAR1/(1d0+VAR*d2wdl2)
   
   !Eq. 18, Update PMU:
   !   self%dPMUdt = self%p%VAR*(self%dmuNetdl - self%dgdlbar - self%dwdl)
   PMU = PMU + log(1D1)  ! Restore to positive values
   if (dmuNetdl .gt. 0d0) then
      PMU1 = PMU  + VAR * dmuNetdl
   else
      PMU1 = PMU/(1d0-VAR/PMU*dmuNetdl)
   endif
   
   PMU1 = PMU1/(1d0 + VAR/PMU*dwdl)
   !Contrain the PMU and VAR: 
   PMU1 = min(max(PMU1,0.01),PMUmax)
   VAR1 = min(max(VAR1,0.01),VARmax)
          
    ! For production/destruction matrix:
      
   pp_ND=dtdays*params(irDN)*DET*tf_z   
   pp_NZ=ZOO*RES        
   pp_DZ=ZOO*EGES+Zmort 
   pp_ZP=ZOO*INGES      
   PPpn =PHY*(muNet + 0.5*VAR*d2muNetdl2)
   pp_PN=max(0.5*( PPpn + abs(PPpn)), 0d0)
   pp_NP=max(0.5*(-PPpn + abs(PPpn)), 0d0)
   
   DET1 = (DET + pp_DZ)/(1d0 + dtdays*params(irDN)*tf_z)

   DET1 = max(DET1,1D-16)
   
   NO31 = (NO3+pp_ND+pp_NZ+pp_NP)/(1d0+pp_PN/NO3)
    
   PHY1 =(PHY+pp_PN)/(1d0+(pp_ZP+pp_NP )/PHY)

   
   ZOO1 = (ZOO+pp_ZP)/(1d0+ EGES+ZOO*dtdays*params(imz)*tf_z+RES)
   
   PMUPHY1 =(PMUPHY+PHY*PMU1+PMU*PHY1)/(1d0+ 2d0*PHY*PMU/PMUPHY)
   VARPHY1 =(VARPHY+PHY*VAR1+VAR*PHY1)/(1d0+ 2d0*PHY*VAR/VARPHY)
   
   Varout(oNO3,k)   = NO31
   Varout(oPHY(1),k)= PHY1
   Varout(oZOO,k)   = ZOO1

   ! Calculate total CHL:
   Varout(oCHLt,k)  = PHY1 / QNAvg * ThetaAvg
   Varout(oPPt,k)   = pp_PN/ QNAvg / dtdays
   Varout(oDET,k)   = DET1
   Varout(oPMU,k)   = PMUPHY1
   Varout(oVAR,k)   = VARPHY1

   if(do_IRON) Varout(oFER,k) = Fe

   Varout(oTheta(1),k) = ThetaAvg
   Varout(oQN(1)   ,k) = QNAvg
   Varout(omuNet(1),k) = PP_PN/dtdays/PHY
   Varout(oGraz(1) ,k) = PP_ZP/dtdays/PHY
   Varout(oZ2N     ,k) = PP_NZ/dtdays
   Varout(oD2N     ,k) = PP_ND/dtdays
   Varout(odmudl   ,k) = dmuNetdl/dtdays
   Varout(od2mu    ,k) = d2muNetdl2/dtdays
   Varout(od2gdl   ,k) = d2gdl2bar/dtdays
   Varout(ow_p(1)  ,k) = w_pAvg  !Should not corrected by dtdays

   ! Calculate size-fractionated Chl
   ! Calculate the proportions of different phyto size fractions
   ! Here the PMU1 is the normal PMU + log(10) (to avoid negative values)
   pCHL(4) = pnorm(PMU_1 +log(1d1),PMU1,sqrt(VAR1))
   pCHL(3) = pnorm(PMU_3 +log(1d1),PMU1,sqrt(VAR1))
   pCHL(2) = pnorm(PMU_10+log(1d1),PMU1,sqrt(VAR1))
   pCHL(1) = 1d0-pCHL(2)
   pCHL(2) = pCHL(2)-pCHL(3)
   pCHL(3) = pCHL(3)-pCHL(4)

   do j = 1, 4
      Varout(oCHLs(j),k) = pCHL(j) * Varout(oCHLt,k) ! Absolute concentration
   enddo

ENDDO
END SUBROUTINE FLEXEFT_CONT 
!========================================================
real function pnorm(x,mean,sd)
! This function calculates the area of the tail left to x 
! of the Gaussian curve with mean and sd
implicit none
real, intent(in) :: x, mean, sd
real             :: x_
real, parameter  :: inv_sqrt2 = 1d0/sqrt(2d0)

x_    = (x-mean)/sd
pnorm = 0.5*(1d0+erf(x_*inv_sqrt2))
return
end function pnorm
!========================================================
real function TEMPBOL(Ea,tC)
implicit none
!DESCRIPTION:
!The temperature dependence of plankton rates are fomulated according to the Arrhenuis equation. 
! tC: in situ temperature
! Tr: reference temperature
!
!INPUT PARAMETERS:
real, intent (in) :: Ea, tC
! boltzman constant constant [ eV /K ]
real, parameter   :: kb = 8.62d-5, Tr = 15D0

TEMPBOL = exp(-(Ea/kb)*(1D0/(273.15 + tC)-1D0/(273.15 + Tr)))
return 
end function TEMPBOL
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
real,    PARAMETER ::EM=-0.367879441171442,&           ! -EXP(-1)
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
       ZN   =ALOG(XX/WAP)-WAP
       TEMP =1D0+WAP
       TEMP2=TEMP+C23*ZN
       TEMP2=2D0*TEMP*TEMP2
       WAP  =WAP*(1D0+(ZN/TEMP)*(TEMP2-ZN)/(TEMP2-2.E0*ZN))
!     END DO
    RETURN
END FUNCTION WAPR

!====================================================
real function ScaleTrait( logsize, star, alpha ) 
implicit none
real, intent(IN) :: logsize, star, alpha

! Calculate the size-scaled value of a trait
! for the given log (natural, base e) of cell volume as pi/6*ESD**3 (micrometers). 

ScaleTrait = star * exp( alpha * logsize )

return
end function ScaleTrait

!====================================================

real function PenAff( logsize, alpha, Pfac, lmin ) 
implicit none
real, intent(IN) :: logsize, alpha, Pfac, lmin 

!A 'penalty' function to reduce the value of affinity for nutrient at very small cell sizes
!in order to avoid modeling unrealistically small cell sizes.  This is needed because affnity
!increases with decreasing cell size, which means that under low-nutrient conditions, without
!such a penalty, unrealistically small cell sizes could be predicted.
!This penalty function becomes zero at logsize = lmin.   
   
  PenAff = 1.0 - exp(Pfac*alpha*(logsize - lmin))
end function PenAff
!====================================================
real function grazing(Hollingtype, Ksat, Prey)
implicit none
real,    intent(in) :: Ksat, Prey
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
    grazing = 1d0-exp(-log(2d0)*Prey/Ksat) 

END SELECT
return
end function
!===================================================
END Module BIO_MOD
