MODULE sub_mod
use Interface_MOD
USE gammaf90
USE mtmod
USE gasdevf90
USE mGf90
implicit none
!!$ declare all variables used within the model itself
public
integer           :: error  = 0, err

! number of samples to be output for the ensemble   
integer           :: EnsLen  = 2 
!Total number of iterations
integer           :: nruns   = 10
! number of runs for initial 'Burn-in'
integer           :: BurnInt = 3 
real              :: alpha13 = 0d0
integer           :: jrun
integer           :: AMacc, DRacc, paramrc, startrun = 0, length, sdwt = 1
! A scratch set of parameters
real, allocatable :: cffpar(:)
real, allocatable :: Apvcurr(:)
real, allocatable :: Apvguess(:)
real, allocatable :: Apvbest(:)
real, allocatable :: Apvmean(:)

!  The arrays of parameter values for the MC chain are defined in terms of only those parameters varied, all are normalized values
real, allocatable :: subpguess(:)
real, allocatable ::   subppro(:)
real, allocatable ::  subppro2(:)
real, allocatable ::  subpbest(:)
real, allocatable ::  subpcurr(:)
real, allocatable ::  subpmean(:)
real, allocatable :: subpcurrmean(:)
real, allocatable ::         sdev(:)
real, allocatable ::    sigmabest(:)                 
real, allocatable ::    sigmamean(:)
real, allocatable ::     BestSSqE(:)
real, allocatable ::       TwtSSE(:), dumE(:)

contains

!-----------------------------------------------------------------------
!Major computation is called here. 
! Calculates model output and Sum of Squares
subroutine CalcYSSE(pars,modval,ssqe)
implicit none

! parameter sets for calculating sum of squares
real, intent(in)  :: pars(NPar)    
real, intent(OUT) :: modval(TNobs)

! Log transformed observational values
real              :: obsval(TNobs)

! Minimal and maximal values for transformation
real              :: max_y,min_y

real, allocatable :: OBS_data(:), MOD_data(:)
! Sum of squares for each data type
real, intent(OUT) :: ssqe(NDTYPE)          

! The model output to match with observational data: 
real              ::  TIN_out(nrow(1),1)
real              ::  CHL_out(nrow(2),1)
real              ::  NPP_out(nrow(3),1)
real, allocatable :: Size_out(:,:)

integer, parameter   :: logtransform = NO
integer, parameter   :: backlog      = YES
integer              :: don,i,k

if (Model_ID == EFTdiscrete .or. Model_ID == EFTcont) then
   allocate(Size_out(nrow(4),4), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
endif

  obsval(:) = 0d0
  modval(:) = 0d0
  
! Make the model read in parameters generated in the main program:      
  params(:) = pars(:)

if (Model_ID == EFTdiscrete .or. Model_ID == EFTcont) then
   call Timestep(TIN_out,CHL_out, NPP_out, Size_out)
else
   call Timestep(TIN_out,CHL_out, NPP_out)
endif
! Match the final year result of the model output to the observational data 

  ssqe(:)   = 0d0
  ! Normalize all the data (including both obs. and mod.) before calculating ssqe 
  !(to give equal weight to each data point)

k = 0
Do i = 1, NDTYPE

  allocate(OBS_data(NDPTS(i)))
  allocate(MOD_data(NDPTS(i)))

  selectcase(i)
  
  case(1)   ! Nitrate
    OBS_data = TINData(:,3)
    MOD_data = TIN_out(:,1)
  case(2)   ! CHL
    OBS_data = CHLData(:,3)
    MOD_data = CHL_out(:,1)
  case(3)   ! NPP
    OBS_data = NPPData(:,3)
    MOD_data = NPP_out(:,1)

  case(4,5,6,7)   ! Size-fractionated Chl
      OBS_data = SizeData(:,i-1)
      MOD_data = Size_out(:,i-3)
  case default
    print *, 'Errors in selecting data types! Quit!'
    stop
  endselect

    max_y =     maxval(OBS_data,1)
    min_y = max(minval(OBS_data,1),0d0) ! Must be positive

! Transform (**0.25) both model and obs. data and normalize between 0 and 1
!subroutine transform(N,x,xmin,xmax,convert,y)
    call transform(NDPTS(i),OBS_data,min_y,max_y,logtransform,  &
         obsval((k+1):(k+NDPTS(i)))  )

    call transform(NDPTS(i),MOD_data,min_y,max_y,logtransform,  &
         modval((k+1):(k+NDPTS(i)))  )

    do don = (k+1),(k+NDPTS(i))

       ssqe(i) = ssqe(i) + (modval(don) - obsval(don))**2 
    enddo
    

    ! Convert back to absolute values
    call transform( NDPTS(i),     modval((k+1):(k+NDPTS(i))),   &
         min_y,max_y,backlog,     modval((k+1):(k+NDPTS(i)))  )
   
    k = k + NDPTS(i)
    
    deallocate(OBS_data)
    deallocate(MOD_data)
  Enddo

end subroutine CalcYSSE
!-----------------------------------------------------------------------
function CalSSQE(Npars) result(SS)
implicit none
! Used to call the model but NOT write output (for most calls in the chain) 
real, intent(IN)  :: Npars(NPar)   ! Normalized parameters

real              :: Ymod(TNobs)
real              :: Apars(NPar)  !Absolute parameters
real              :: SS(NDTYPE)
  ! Convert normalized parameters to real parameters:
  Apars = Apv_(Npars)

  ! Calculate SSqE based on parameter input (a real model run):
  call CalcYSSE(Apars,Ymod,SS)

end function CalSSQE
!-----------------------------------------------------------------------
! Used to call the model with best parameters AND write output to the "bestout" file
subroutine model(fint, subp, SS)
implicit none

! Unit of Data file to be written
integer, intent(IN)  :: fint

! Current parameters
real,    intent(IN)  :: subp(Np2Vary)

! Calculated SSQE
real,    intent(OUT) :: SS(NDTYPE)

real    :: Ymod(TNobs) 
integer :: don

! assign Apv at each step
  Apv=Apv_(subp)

  call CalcYSSE(Apv,Ymod,SS)
  write(fint, 3000) '  DOY    ',    &
                    'Depth    ',    &
                    'Data_type',    &
                    'Mod_Value'
   ! Write the best simulation results into the file 
  do don = 1, TNobs
     write(fint, 3100) OBS_DOY(don),OBS_Depth(don),OBS_Label(don),Ymod(don)
  enddo
      
3000   format(1x,     3(A10), 3x, A10)
3100   format(1x, 2(F6.1, 2x), 5x, A5, 1x, 1pe12.3 )
end subroutine model
!-----------------------------------------------------------------------
subroutine modelensout(efint, runnum, subp, SS)
  ! Used to call the model with the current parameter values 
  ! AND write output to the file for the Ensemble of model output 
  ! (only once every "outputint" runs)  
integer, intent(IN)  :: efint, runnum
real(8), intent(IN)  :: subp(Np2Vary)
real(8), intent(OUT) :: SS(NDTYPE)
real(8)              :: Ymod(TNobs)
integer              :: don

  Apv = Apv_(subp)
  call CalcYSSE(Apv,Ymod,SS)

  do don = 1, TNobs
     write(efint, 3100) runnum, OBS_DOY(don),OBS_Depth(don),&
                                OBS_Label(don), Ymod(don)
  enddo

3100   format(1x, I10,1x,  2(F6.1, 2x), A10, 1x, 1pe12.4)
end subroutine modelensout

!-----------------------------------------------------------------------
subroutine modelnooutput(subp, SS)
implicit none
! Used to call the model but NOT write output (for most calls in the chain) 
real, intent(IN)  :: subp(Np2Vary)
real, intent(OUT) :: SS(NDTYPE)
real              :: Ymod(TNobs)
  ! The subp is the normalized value! 
  Apv=Apv_(subp)
  call CalcYSSE(Apv,Ymod,SS)

end subroutine modelnooutput

!-----------------------------------------------------------------------
subroutine write_bestsigma
implicit none
integer :: i,k

! Write into best sigma file:
open(bsfint, file=bsfn, action='write',status='replace' )
write(bsfint,1200) 'At Run #        ',jrun
write(bsfint,1210) 1D2*real(AMacc)/real(jrun-startrun),       &
                   1D2*real(DRacc)/real(jrun-startrun)
write(bsfint,1220) 'Best LogL =    ', BestLogLike

! To calculate weighted SSqE(SS/sigma**2, Eq. 20), just to compare to previous non-adaptive MCMC runs
TwtSSE = 0
do k = 1, NDTYPE
   TwtSSE = TwtSSE + BestSSqE(k)/(sigmabest(k)**2)
enddo
write(bsfint,1220) 'Best weighted SSqE = ', TwtSSE
write(bsfint,1300)  NewLogLike, CurrLogLike, BestLogLike
write(bsfint,1400) (SigmaLabel(i), sigmabest(i), sigmamean(i), i = 1, NDTYPE)
close(bsfint)
1400 format('                  sigmabest              sigmamean ',/, &
           100(a15,2(1x,1pe20.3),/) )
1200 format(a15,1x,i16)
1210 format('** % 1st Accept. = ',1x,1f8.2,                         &
            '     ** % 2nd Accept. = ',1f8.2)

1220 format(a25,1x,1pe11.3)
1300 format('*** LogL:    New ',1pe11.3,'     Curr ',1pe11.3,  &
           '     Best ',1pe11.3)
1330 format('***        with  Subpcurr                Subpguess ',   &
           '    -->       Subpbest ','               Subpmean ',/,   &
           100(a15,4(1x,1pe20.3),/) )

end subroutine 
!-----------------------------------------------------------------------
subroutine write_bestpar
implicit none
integer :: i
! Write into best parameter file:
open(bpfint, file=bpfn, action='write',status='replace'  )
write(bpfint,1200) 'At Run #        ',jrun
write(bpfint,1210) 1D2*real(AMacc)/real(jrun-startrun),         &
                   1D2*real(DRacc)/real(jrun-startrun)

write(bpfint,1220) 'Best logL =    ', BestLogLike
write(bpfint,1300)  NewLogLike, CurrLogLike, BestLogLike

Apvcurr = Apv_(subpcurr)
Apvguess= Apv_(subpguess)
Apvbest = Apv_(subpbest)
Apvmean = Apv_(subpmean)

write(bpfint,1330) (ParamLabel(i),   &
  Apvcurr(i),                        &
  Apvguess(i),                       &
  Apvbest(i),                        &
  Apvmean(i), i = 1, Np2Vary)
close(bpfint)

1200 format(a15,1x,i16)
1210 format('** % 1st Accept. = ',1x,1f8.2,                         &
            '     ** % 2nd Accept. = ',1f8.2)

1220 format(a25,1x,1pe11.3)
1300 format('*** LogL:    New ',1pe11.3,'     Curr ',1pe11.3,  &
           '     Best ',1pe11.3)
1330 format('***        with  Subpcurr                Subpguess ',   &
           '    -->       Subpbest ','               Subpmean ',/,   &
           100(a15,4(1x,1pe20.3),/) )

end subroutine write_bestpar
!-----------------------------------------------------------------------
function Invmat(C) result(InvC)
implicit none
real, intent(in)  ::  C(NPar*(NPar+1)/2)
real              ::  InvC(NPar,NPar)

real              ::  work(NPar), U(NPar*(NPar+1)/2)
!!$ Invert the symmetric covariance matrix 
  call syminv(C, NPar, U, work, nullty, error )

!!$ Unpack the compressed inverse matrix U into the full (2-D) inverse covariance matrix
  call unpack(NPar,U,InvC)
end function
!-----------------------------------------------------------------------
function newsigma(SS) result(sigma_)
implicit none
integer             :: i
real, intent(in)    :: SS(NDTYPE)
real                :: sigma_(NDTYPE)
!!$  sigma is randomly generated every time based on CurrSSqE
    ! Eq. 22 in Laine 2008
    do i = 1, NDTYPE
       sigma_(i) = sqrt( 1d0/  &
       gammadev( (n0d2(i) +NDPTS(i)/2d0),   &
                 1d0/(n0d2(i)*S02(i) + SS(i)/2d0) ) )
    enddo
end function newsigma

!-----------------------------------------------------------------------
logical function check_bounds(pars) result(NOTOK)
implicit none
real,    intent(in)  :: pars(NPar)
integer              :: k

NOTOK = .FALSE.
do k = 1, NPar  ! Loop over the Parameters to be varied

   if ((pars(k) .lt. MinValue(k)) .or. (pars(k).gt.MaxValue(k))) then
      
      NOTOK = .true.
      exit
   endif
        
enddo ! loop over Parameters 
return
end function

!-----------------------------------------------------------------------
function goodPar(Rchol,Pars) result(y)
implicit none
! Rchol is the lower triangluar  Cholesky factor of
! the Proposal Covariance Matrix, 
! which is calculated only every outputint steps 
! (it remains constant for each set of outputint steps). 
real,  intent(in)   :: Rchol(NPar*(NPar+1)/2)

! The desired mean of the new params
real,  intent(in)   :: Pars(NPar)
real  :: y(NPar)

! Absolute parameters
logical             :: ISBADPAR

! Pars is a one dimensional array holding proposed values for the subset of parameters being varied. 
    ISBADPAR = .true.

! If any parameter gets outside it's range, start over to generate a 
! different parameter set

    DO WHILE (ISBADPAR) 
    ! generate new set of parameters
       y = multiGauss(Rchol,Pars,NPar)
      
! Restrict all the parameters within the bounds: 
       cffpar  = Apv_(y)

       ISBADPAR = check_bounds(cffpar)
    ENDDO
End function
!-----------------------------------------------------------------------
function NewPAR(Pcvm_, OldPar) result(y)
implicit none
real, intent(in)   :: OldPar(NPar)
real, intent(in)   ::  Pcvm_(NPar*(NPar+1)/2)

real  :: y(NPar)

! Calculate the Cholesky factor for Pcvm, which is called Rchol here. 
call cholesky(Pcvm_,NPar,NPar*(NPar+1)/2,  &
                    Rchol,nullty,error)

y = goodPar(Rchol, OldPar)

end function
!-----------------------------------------------------------------------
subroutine write_ensemble
implicit none
integer i

cffpar = Apv_(subpcurr)

! Write Ensemble file(s) (Run #, LogLike and Parameters)
write(epfint,1850) jrun, CurrLogLike, (cffpar(i), i = 1, Np2Vary)
write(esfint,1850) jrun, CurrLogLike, &
     (sigma(i),  i = 1, NDTYPE),      &
     (CurrSSqE(i),i= 1, NDTYPE)

1850 format(i9,1x,100(1pe12.3,2x))
end subroutine
!-----------------------------------------------------------------------
!!!
!!! second stage DR acceptance probability
!!! assumes Gaussian proposals, global variable iC = inv(C) is the inverse of
!!! first stage proposal of the point that was rejected
!!! 
subroutine MCMC_DR_alpha13(C, &
                           par1,      logLike1, &
                           par2,      logLike2, &
                           par3, ss3, logLike3, alpha13)

implicit none
real, intent(in) :: C(NPar*(NPar+1)/2)      ! the first stage proposal

! The SSqE of  new2 params
real, intent(in) ::  ss3(NDTYPE)  

! Three sets of parameters
real, intent(in) :: par1(NPar), par2(NPar), par3(NPar)

! The probability of accepting the first move (already rejected)
real, intent(in) :: logLike1, logLike2

! The probability of the second move given the current position 
! and the first move (already rejected)
! logLike3 is the loglikelihood of the second move
real, intent(out):: alpha13, logLike3

real             :: alpha12
real             :: l2, q1, alpha32
real             :: iC(NPar, NPar)         ! the uncompacted form of inverse of C  
real             :: cff1(NPar), cff2(NPar),dpar32(NPar),dpar12(NPar)
!integer          :: m,q
! Calculate iC:
iC      = Invmat(C)

!  write(6,*) ' The Inverse of the Covariance Matrix: '
!
!  do m = 1, NPar
!     write(6,1100) (IC(m,q), q = 1, NPar)
!  end do
!
!1100 format(<NPar>(1pe9.1,1x))
! The probability of accepting the first move
alpha12 = min(1d0, exp(logLike2-logLike1))

! Calculate the logLike of the second move
     
    logLike3 = CalcLogLike(ss3, sigma, par3)

! Calculate the probability (alpha32) of accepting y2 given y3

  if (alpha12 == 0d0) then
     alpha32 = 0d0
  else

! the loglikelihood from y3 to y2
!! oldpar is par3 and oldSS is ss3!!!

     alpha32 = min(1d0, exp(logLike2 - logLike3))
  endif

  ! l2 = log(posterior(par3) - posterior(par1)  )
  
  l2 = logLike3-logLike1

  dpar32 = par3-par2
  dpar12 = par1-par2

  call matmuls(iC, dpar32, cff1) 
  call matmuls(iC, dpar12, cff2) 

!  q1 = -0.5d0*( &
!       sum(  matmuls(iC,(par3-par2)) * (par3-par2) ) - &
!       sum(  matmuls(iC,(par1-par2)) * (par1-par2) ) )

 q1 = -0.5d0*(sum(cff1*dpar32) - sum(cff2 * dpar12))


  ! alpha2(x,y1,y2)
  alpha13 = min(1d0, exp(l2+q1)*(1d0 - alpha32)/(1d0 - alpha12))

end subroutine MCMC_DR_alpha13
!-----------------------------------------------------------------------
subroutine MCMC_adapt(BURNIN)
! Major MCMC subroutine
implicit none
logical, intent(in)       ::  BURNIN

real   :: SS3(NDTYPE), alpha13
real   :: NewLogLike2
real   :: PCVM2(NPar*(NPar+1)/2)
! Output also includes:
! 1) Params, CVM, PCVM, SSqE, sigma, loglikehood at each step
! 2) Best of the above quantities so far
integer :: i,j, k, EnsLen1,intv,jrun1

! Downscaling factor for Delayed Rejection MCMC
real,    parameter        :: DRScale = 1d-2

! The number of times of adaptation (i.e. number of ensembles)

If (BURNIN) Then
   EnsLen1   = 1
   intv      = BurnInt
Else
   EnsLen1   = EnsLen
   intv      = (nruns-startrun)/EnsLen
Endif
! Acceptance counts for Adaptive MCMC
AMacc = 0

! Acceptance counts for Delayed Rejection MCMC
DRacc = 0

BestLogLike = -1d12  ! A very large, negative number for very low probability
MeanLogLike = 0d0

DO i = 1, EnsLen1

 if (.not. BURNIN) then
!  write output to Ensemble files for simulated values
   call modelensout(eofint,jrun,subpcurr, dumE)

! Write status file(s) (snapshot of: Run #, Cost, acceptance and Parameters)
! Open & Close each time, so that it is written even if the program crashes
   call write_bestpar
   call write_bestsigma

!  write best results into the best file
   open(bofint, file=bofn, action='write',status='replace')
   savefile=.FALSE.
   call model(bofint, subpbest, dumE ) !dumE: SSqE for best parameters
   close(bofint)
 endif
! end of block of outputting code
   subcycle: DO j = 1, intv
     
      if (.not. BURNIN) then
          jrun1 = jrun -startrun
       subpmean =((real(jrun1)-1d0)*subpmean + subpcurr)    /real(jrun1)
       sigmamean=((real(jrun1)-1d0)*sigmamean+    sigma)    /real(jrun1)
     MeanLogLike=((real(jrun1)-1d0)*MeanLogLike+CurrLogLike)/real(jrun1)
   ! write current parameters and sigma into output file
      Call write_ensemble

      endif
         call modelnooutput(subppro, SSqE)

      !   New log-likelihood
      NewLogLike = CalcLogLike(SSqE,     sigma, subppro)

      !   Old log-likelihood
      CurrLogLike = CalcLogLike(CurrSSqE, sigma,subpcurr)

      IF( log(grnd())  .lt. (NewLogLike - CurrLogLike) ) THEN
      !acc is a tally of the acceptance rate, the rest are parameters
           AMacc       = AMacc+1
           CurrLogLike = NewLogLike
           CurrSSqE    = SSqE
           subpcurr    = subppro
      ELSE
      ! DR:
      ! Propose a second move (Y2) based on scaled proposal covariance matrix
      ! and current position
           Pcvm2      = DRScale*Pcvm
           subppro2   = NewPAR(Pcvm2, subpcurr)
         !      Calcuate the new SSqE
           SS3        = CalSSQE(subppro2)

         ! Calculate the acceptance probability of the second move
           call MCMC_DR_alpha13(Pcvm,      &
                  subpcurr,  CurrLogLike,  &
                  subppro,    NewLogLike,  &
                  subppro2,      SS3,  NewLogLike2, alpha13)

           ! Judge whether the second move should be accepted or not
           if( GRND() .lt. alpha13 ) then
               !accept the second move (Y1)
              DRacc       = DRacc+1
              CurrLogLike = NewLogLike2
              CurrSSqE    = SS3
              subpcurr    = subppro2
           endif
      ENDIF

      !	this is to keep track of the 'best' run which 
      !	is not strictly speaking the purpose of the
      !	assimilation but an interesting output
      !
      IF( NewLogLike .GT. BestLogLike ) THEN
         BestLogLike = NewLogLike     
         BestSSqE    = SSqE
         sigmabest   = sigma
         subpbest    = subppro
      ENDIF

! Update the Covariance matrix for the Parameter values, using the recursion formula
! assuming a weight of 1 for each step, updating at every step. 
        
      call UpdateCVM(Cvm,subpcurrmean, Rwtold, subpcurr, subpcurrmean, Cvm)

      ! increment the wt. for the existing covariance matrix
      Rwtold = Rwtold + 1d0 !Weight of Cvm (Covariance matrix)

      !   sdwt   = sdwt   + 1 

      ! Propose newpar based on old par
       subppro = NewPAR(Pcvm, subpcurr)

      ! Update sigma:
       sigma = newsigma(CurrSSqE)

      ! Update the counter for number of simulations
      jrun= jrun+ 1

   ENDDO subcycle ! END of loop of each cycle 

! Calculate new PCVM based on CVM
! Set the Proposal covariance matrix, by scaling the Parameter Covariance matrix
   Pcvm = Cvm*Spcvm/NPar
! Add a small term to the diagonal entries, so that the matrix will not be singular. 
   do k = 1, NPar
     Pcvm(k*(k+1)/2)=Pcvm(k*(k+1)/2) + CvEpsilon*Spcvm/NPar
   enddo
ENDDO

End subroutine MCMC_adapt
!-----------------------------------------------------------------------
END MODULE sub_mod
