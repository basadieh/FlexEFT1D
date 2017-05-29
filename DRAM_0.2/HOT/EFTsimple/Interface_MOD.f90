MODULE Interface_MOD
! User defined section:
! The read data section better to be put together with the modeling section
! because the model needs to interpolate the model output to match with the OBS. data
USE MOD_1D
USE syminvf90
IMPLICIT NONE
public
! A scratch set of parameters
real, allocatable :: cffpar(:)
real, allocatable :: Apvcurr(:)
real, allocatable :: Apvguess(:)
real, allocatable :: Apvbest(:)
real, allocatable :: Apvmean(:)

!  The arrays of parameter values for the MC chain are defined in terms of only those parameters varied, all are normalized values
real, allocatable ::    subpguess(:)
real, allocatable ::      subppro(:)
real, allocatable ::     subppro2(:)
real, allocatable ::     subpbest(:)
real, allocatable ::     subpcurr(:)
real, allocatable ::     subpmean(:)
real, allocatable :: subpcurrmean(:)
real, allocatable ::         sdev(:)
real, allocatable ::    sigmabest(:)                 
real, allocatable ::    sigmamean(:)
real, allocatable ::     BestSSqE(:)
real, allocatable ::       TwtSSE(:)
real, allocatable ::         dumE(:)

! Weights and Sums of Squares 
! For the standard deviations (sigma) used to give weights to each type of data. 
character(LEN=10), allocatable :: SigmaLabel(:)  ! to be assigned by the subroutine SetSigmaLabels
character(LEN=10), allocatable :: SSqELabel(:)

!A standard error (sigma) of the data will be estimated by Gibbs sampling for each data type, respectively. 
! sigma declared here.
real, allocatable :: sigma(:)
! Sums of square errors, for each data type, respectively
real, allocatable :: SSqE(:)
real, allocatable :: CurrSSqE(:)
real, allocatable :: SS3(:) 
      
! Parameters for the prior component of the SE, assumed to have a Gamma distribution
! *** Note: these should be set differently, depending on the data set. 
! *** Ideally, they should represent the distribution of MINIMUM observation 
!     and sampling error (and perhaps also minimum model-data mismatch). 
! *** Using Gibbs sampling, the algorithm will sample the distribution of each sigma, 
!     which will include the contributions of BOTH model-data mismatch AND actual
!     observation and sampling error.
! accuracy of the prior for sigma (# of observations averaged for each data pt.) 

! ideally, the number of replicates for each pt. 
real, allocatable :: n0d2(:)

! the mean of the prior for sigma (standard error squared) 
real, allocatable :: S02(:)  ! lower values give more weight to data

!!$ Model-specific parameters 
      
! To create a vector containing only the subset of parameters to be varied 
! in the Monte Carlo chain (this can be any subset of the NPar parameters) 
! The number of Parameters to be varied (fitted), including VAR 
integer              :: Np2Vary      

! Parameter Values, the quantities to be varied by the MCMC algorithm. 
! *** Note: the subroutines in submodel.f90 should be written to access the local parameter
! values as, e.g., the current value of parameter "theta1" = apv(pitheta1)
! This allows the main program to pass different parameter sets to the subroutines, 
! for example the newly proposed parameter set, or the best-fit parameter set. 

real, allocatable :: Apv(:)   ! Actual parameter values 
real, allocatable :: Npv(:)   ! Normalized parameter values [dimensionless, for AM algorithm]
real, allocatable :: Ipv(:)   ! Initial parameter values

!!$ Parameter Scaling factor for normalizing (non-dimensionalizing) equations 
! Fixed throughout the program?
!real, allocatable :: pSF(:)
      
!  Bounds on the ranges of parameter values are defined using the master index (of pointers)
real, allocatable :: MaxValue(:), MinValue(:)

!!$ Variables for the distributions of parameter values, priors, etc. 

! Priors (prior estimates) for parameters: 
! Assume a Gaussian distribution, independenctly for each parameter, respectively. 
! In this code, these are set in a subroutine (below). 
real, allocatable :: nuPar(:)   ! mean of each parameters prior
real, allocatable :: etaPar(:)  ! standard deviation of each parameters prior distribution. 

!For estimating priors, assumed CV for all prior parameter estimates
!Note: Lower values of CVprior will give more weights to the Prior Estimates 
!in the Log Likelihood; i.e. the fit will be penalized more (lower LogL) for 
!parameter values deviating from those estimates. 
!**** CVprior specifies the width (coefficient of variation) for the prior distribution of each parameter
!**** Too wide a distribution will result in extremely low acceptance rates in the AM algorithm, 
!**** especially if parameter sets are rejected whenever ANY proposed param. value is < 0. 
!**** Consider a random value x, chosen from  Gaussian distribution of mean 1 and specified CV:
!****   CV       P[ x > 0 ]   P[ x > 0 ]^10
!****   0.25       0.99997       0.9997 
!****   0.333      0.99866       0.987 
!****   0.5        0.97725       0.794
!***    1.0        0.84134       0.178 
!*** The third column is the probability that 10 independently drawn x will all be > 0.  
!****  As the number of parameters being fitted increases, 
!      the chances of rejecting a param set increase, which will require smaller CVprior. 
real, parameter   :: CVprior = 2d0    

! Covariance matrix for Priors, initialized to zero
real, allocatable :: PriorCvm(:)

! Inverse Prior Covariance matrix (NOT compacted!)  
real, allocatable :: InvPriorCvm(:,:)

! Logs of posterior probabilities with New (proposed) and current values of parameters
real  :: BestLogLike, CurrLogLike, MeanLogLike, NewLogLike, NewLogLike2

! Cholesky factor for the matrix of covariances of the parameters
! (lower triangular matrix, stored by row, as a one dimensional array) 
real, allocatable :: Rchol(:)

! Covariance matrix, initialized to zero
real, allocatable :: Cvm(:)

! Proposal Covariance matrix, initialized to zero 
real, allocatable :: Pcvm(:)

! For DRAM
real, allocatable :: Pcvm2(:)

! for storing the square root of the determinant of Pcvm 
real              :: rtDetPcvm                  

! A scaling parameter to make certain that the proposal covariance matrix does not become singular: 
real, parameter   :: Spcvm = 2.4d0*2.4d0

! A small term to be added to the diagonal 
real, parameter   :: CvEpsilon = 1d-16

! weight for initial guess of proposal covariances  
real              :: Rwtold = 1d0                

!!$ Gibbs sampling (calculates weighs for data) 

! Marko Laine (PhD thesis, Lapeenranta University of Technology, Finland, 1998) 
! writes the parameters for the distribution as Gamma( n0d2, S02*n0d2 ), 
! such that the mean of the distribution for sigma = S02. 
! (With the parameters written that way, mean( Gamma(k, theta) ) = k / theta).  
! He assumes that sigma**-2  is distribtuted as 
! Gamma( k = (n0d2 + n/2), theta =(S02*n0d2 + SSn/2) ). 
! **Be careful when putting parameters into the function to generate Gamma distributed 
! pseudo-random numbers,because that function, gamamdev, is written to take parameters 
! alpha and beta such that:
!                mean( Gamma(Alpha, Beta) ) = Alpha * Beta. 
! (Both ways of writing the Gamma distribution are common.) 
! Thus, using the routine gammadev to generate the random samples from the distribution 
! for sigma, parameters should input as:
!          gammadev( Alpha = (n0d2 + n/2), Beta =( 1/(n0d2*S02 + SSn/2) )
! that is, the second parameter is input as the reciprocal of what Laine (2008) wrote. 

!!$ For Input/Output (file names, indeces, etc). 
  ! indeces and labels
  ! for variables used to interface the main program with the Phytoplankton model
integer, parameter :: Yes = 1, No = 0
integer            :: First = Yes
integer            :: nullty  ! i/o for cholesky subroutine

!!$ Output files

! best file for model output
integer, parameter :: bofint = 25

! best file for parameter values (status file)
integer, parameter :: bpfint = 26 

! best file for sigma values 
integer, parameter :: bsfint = 27 

! file for reading best parameter values
integer, parameter :: rpfint = 28 

!ensemble file for model output
integer, parameter :: eofint = 29

! ensemble file  for parameter values
integer, parameter :: epfint = 30 

! ensemble file  for sigma values 
integer, parameter :: esfint = 31

!!$ Output filenames
character(7), parameter :: bofn = 'bestout'

! Output files for best parameters
character(8),  parameter :: bpfn = 'bestpar'
character(8),  parameter :: bsfn = 'bestsig'  !One file
character(8),  parameter :: rpfn = 'Best_par' !One file
character(8),  parameter :: eofn = 'ensout'
character(20), parameter :: epfn = 'enspar'
character(20), parameter :: esfn = 'enssig'   !One file

!!$ Mathematical constants for calculations
real, parameter         :: rt2Pi = sqrt(2D0*pi) 

CONTAINS
      
! A function that converts integer to string
character(len=20) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine SetUpArrays  ! sets up the arrays of models and parameters 
implicit none
integer  :: k
!!$  ***  Every Parameter, including those not varied, should be assigned a label 
!!$  ***  CONSISTENT with the indeces for all parameters (piXXX) defined above  

! Select biological model
call choose_model

! Read observational data
call Setup_OBSdata

Np2Vary = NPar

! Initialize all the relevant variables
allocate(Apv(NPar))
allocate(Npv(NPar))
allocate(Ipv(NPar))
allocate(MaxValue(NPar))
allocate(MinValue(NPar))
allocate(nuPar(Np2Vary))
allocate(etaPar(Np2Vary))

! Covariance matrix for Priors
allocate( PriorCvm(Np2Vary*(Np2Vary+1)/2) )
PriorCvm(:) = 1d-16

! Inverse Prior Covariance matrix (NOT compacted!)  
allocate( InvPriorCvm(Np2Vary,Np2Vary) )
InvPriorCvm(:,:) = 1d0

! Cholesky factor for the matrix of covariances of the parameters
! (lower triangular matrix, stored by row, as a one dimensional array) 
allocate(  Rchol(Np2Vary*(Np2Vary+1)/2) )
Rchol(:) = 0d0

! Covariance matrix, initialized to zero
allocate( Cvm(Np2Vary*(Np2Vary+1)/2) )
Cvm(:)   = 0d0

! Proposal Covariance matrix, initialized to zero 
allocate( Pcvm(Np2Vary*(Np2Vary+1)/2) )
Pcvm(:)  = 0d0

allocate(Pcvm2(Np2Vary*(Np2Vary+1)/2) )
Pcvm2(:) = 0d0

allocate(subpguess(Np2Vary))
subpguess(:) = 1d0
allocate(  subppro(Np2Vary))
subppro=subpguess
allocate( subppro2(Np2Vary))
subppro2=subpguess
allocate( subpbest(Np2Vary))  
subpbest=subpguess
allocate( subpcurr(Np2Vary))  
subpcurr=subpguess
allocate( subpmean(Np2Vary))  
subpmean=subpguess
allocate(   cffpar(Np2Vary))  
cffpar=subpguess
allocate( Apvguess(Np2Vary))  
allocate(  Apvbest(Np2Vary))  
allocate(  Apvcurr(Np2Vary))  
allocate(  Apvmean(Np2Vary))  
allocate(subpcurrmean(Np2Vary))  
subpcurrmean=subpguess
allocate(sdev(Np2Vary*(Np2Vary+1)/2)) 

allocate(   sigma(NDTYPE))
sigma(:) = 1d-2
allocate(    SSqE(NDTYPE))
SSqE(:)  = 0d0
allocate(CurrSSqE(NDTYPE))
CurrSSqE(:) = 0d0
allocate(     SS3(NDTYPE))
SS3(:) = 0d0
allocate(    n0d2(NDTYPE))
n0d2(:)= 1d0
allocate(     S02(NDTYPE))
S02(:) = 1d-1
allocate(SigmaLabel(NDTYPE))
SigmaLabel(:)= 'SigmaLabel'
allocate(SSqELabel(NDTYPE))
SSqELabel(:)='SSqELabel'
allocate(sigmabest(NDTYPE))   
sigmabest = sigma
allocate(sigmamean(NDTYPE))   
sigmamean = sigma
allocate( BestSSqE(NDTYPE))    
BestSSqE=SSqE
allocate(   TwtSSE(NDTYPE))   
TwtSSE=SSqE
allocate(     dumE(NDTYPE))
dumE=SSqE

!===================================================
! Maximum and minimum allowed values for each parameter set here  
! using the arrays MaxValue and MinValue. 
!===================================================

! Common parameters:

! Based on Edwards et al. (2015) who gave aI0 ranging from
! 0.001 to 1 (uE m-2 s-1)-1 d-1, corresponding to the values times *4.6/(50/12)
! to get the unit: (mol C gChl-1 (W m-2)-1 d-1)

! Fasham (DSR1 1995) gave 0.025~0.14 unit: (W m-2)-1 d-1
MaxValue(iaI0     ) =  1.1       ! 1*4.6/(50/12)
MinValue(iaI0     ) =  1.1d-3

! Fennel et al. (2006): 0.56~3.5
! Chai et al. (2002): 0.25~0.5
! Li et al. (2001): 0.003
! Chen et al. (2014) gave 0.02~7.4 ug Chl/L, 
! approximating 0.01~4.7 uM nitrogen
! Hansen et al. (1997) gave 0.4~40 ug Chl/L, approximating 0.2 ~ 20 uM N
MaxValue(ikp      ) =  4.7
MinValue(ikp      ) =  0.01

! Fennel et al. (2006): 0.01~0.25
! Kishi et al. (2007): 0.05~0.2
MaxValue(irdn     ) =  0.3
MinValue(irdn     ) =  0.01

! Fennel et al. (2006) gave range of 0.009-25 m/d
! Kishi et al. (2007) gave sinking rate of POC of 40 m/d
MaxValue(iwDET    ) =  25.0
MinValue(iwDET    ) =  0.009

! Fennel et al. (2006) gave 0.025, but gave a range of 0.05~0.25
! Lima and Doney (2004) gave 0.25
! Chai et al. (2002) gave 0.05 for mesozooplankton
! So I decide to give a relatively large range
MaxValue(imz      ) =  0.3
MinValue(imz      ) =  0.02

! Fennel et al. (2006) gave a range of 0.5~1.0, but used 0.6
! Lima and Doney (2004) gave 2.75 to 3.75
! Hansen et al. (1997) gave a range of 10**(-2) ~ 10**(0) (unit: h-1)
! after conversion to d-1, 0.2 ~ 2.4
MaxValue(igmax    ) =  3d0
MinValue(igmax    ) =  0.1

! Model-specific parameters:
select case(model_ID)

case(Geidersimple)
  !Based on the lab dataset from Chen and Laws (in revision)
  !Growth rate normalized to 15 ºC based on linear regression
  !0.025% and 0.975% quantiles
  MaxValue(imu0     ) =  2.35
  MinValue(imu0     ) =  1d-2

!  MaxValue(iPHYini  ) =  5d0
!  MinValue(iPHYini  ) =  1d-8

!  MaxValue(iZOOini  ) =  5d0
!  MinValue(iZOOini  ) =  1d-8
!
!  MaxValue(iDETini  ) =  5d0
!  MinValue(iDETini  ) =  1d-8


  ! This is redfield N:C ratio, bounds based on Geider et al. (2002)
  ! and our own data on Platymonas
  MaxValue(iQ0N     ) =  0.3
  MinValue(iQ0N     ) =  0.05

  ! Fennel et al. (2006): 0.007~1.5
  ! Chai et al. (2002): 0.05~1
  ! Franks (2009): 0.005~3
  MaxValue(iKN      ) =  3d0
  MinValue(iKN      ) =  0.005


!-------------------------------
case(EFTsimple)
!  Following Pahlow et al. (2013), fix mu0 and V0N as 5 per day
!  Net growth rate regulated by A0N and Q0N

!  MaxValue(imu0     ) =  8d0
!  MinValue(imu0     ) =  0.1d0

!  MaxValue(iPHYini  ) =  5d0
!  MinValue(iPHYini  ) =  1d-4
!
!  MaxValue(iZOOini  ) =  5d0
!  MinValue(iZOOini  ) =  1d-4
!
!  MaxValue(iDETini  ) =  5d0
!  MinValue(iDETini  ) =  1d-4

  ! This Q0N is the minimal QN value (unit: mol N: mol C), 
  ! different from Geidersimple model
  ! Litchman et al. (2007) gave a range from 0.01 to 0.07
  ! Pahlow et al. (2013) gave a range from 0.05 to 0.13
  MaxValue(iQ0N     ) =  0.15
  MinValue(iQ0N     ) =  0.01

  ! Maximal rates of N uptake
  !MaxValue(iV0N     ) =  2d1
  !MinValue(iV0N     ) =  1d-8

  select case(nutrient_uptake)
  case(1)
    MaxValue(iKN    ) =  3D0
    MinValue(iKN    ) =  0.005
  case(2)
    ! Pahlow et al. (2013) gave a range of 60 ~ 1000 (unit: m3 molC-1 d-1)
    ! Smith et al. (2015) estimated as 0.15 (unit: m3 mmolC-1 d-1)
    MaxValue(iA0N   ) =  1d0
    MinValue(iA0N   ) =  0.05
  case default
    write(6,*) 'Nutrient uptake option incorrect! Quit!'
    stop
  end select

!-------------------------------
case(EFTdiscrete,EFTcont)  

  if (model_ID .eq. EFTcont) then
     MaxValue(ithetamin) = 0.1
     MinValue(ithetamin) = 0.01
     MaxValue(iQNmin)    = 0.05
     MinValue(iQNmin)    = 0.005
     MaxValue(iPenfac)   = 1D1
     MinValue(iPenfac)   = 1D-1
     
  endif
  MaxValue(ialphaI  ) =  0.1
  MinValue(ialphaI  ) =  -0.5

  select case(nutrient_uptake)
  case(1)
    MaxValue(iKN)     =  3D0
    MinValue(iKN)     =  0.005
    MaxValue(ialphaK) =  0.1
    MinValue(ialphaK) =  -0.5
  case(2)
    MaxValue(iA0N)    =  1d0
    MinValue(iA0N)    =  0.05
    MaxValue(ialphaA) =  0.1
    MinValue(ialphaA) =  -0.5
  case default
    write(6,*) 'Nutrient uptake option incorrect! Quit!'
    stop
  end select
  !MaxValue(iV0N     ) =  8d0
  !MaxValue(ialphaV  ) =  0.4
  MaxValue(iQ0N     ) =  0.15
  MinValue(iQ0N     ) =  0.01

  MaxValue(ialphaQ  ) =  -(1d-9)
  MinValue(ialphaQ  ) =  -0.4

  MaxValue(ialphaG  ) =  2d0
  MinValue(ialphaG  ) =  1d0

  !MinValue(imu0     ) =  2d0
  !MinValue(iV0N     ) =  2d0
  !MinValue(ialphaV  ) =  -0.4

  !MaxValue(iPHYini  ) =  1d2
  !MinValue(iPHYini  ) =  1d-8
  !MaxValue(iZOOini  ) =  1d2
  !MinValue(iZOOini  ) =  1d-8
  !MaxValue(iDETini  ) =  1d2
  !MinValue(iDETini  ) =  1d-8

case default
  write(6,*) 'Model option incorrect! Quit!'
  stop
end select  

! Set the prior parameters
do k = 1, NPar
     Apv(k)   = (MaxValue(k)+MinValue(k))/2d0
     Ipv(k)   = Apv(k)
     Params(k)= Apv(k)
enddo

! Normalize the actual parameter values (Apv), 

!It is essential to normalize the parameter values in order to calculate LogLikelihoods 
!that are meaningful even in cases where parameter values differ greatly within the 
!set of parameters (i.e., p1 is very small, and p2 is very large, in absolute magnitude). 
!In such cases, using the absolute parameter values may give too much weight to those 
!parameters having large absolute values, in the calculation of the prior 
!contribution to the LogLikelihood. 

  Npv = Npv_(Apv)
END SUBROUTINE SetUpArrays
!-----------------------------------------------------------------------
! This function normalize the absolute paramter values
function Npv_(Apv_) result(vals)
implicit none
real, intent(in ) :: Apv_(NPar)
real :: vals(NPar)
integer           :: i

   do i = 1, NPar 
     ! vals(i) = (Apv_(i) - MinValue(i))/(MaxValue(i) - MinValue(i))
      vals(i) = Apv_(i)/Ipv(i)
   enddo

end function Npv_
!-----------------------------------------------------------------------
! This function back calculates the normalized parameter values to absolute paramter values
function Apv_(Npv_) result(vals)
implicit none
real, intent(in ) :: Npv_(NPar)
real              :: vals(NPar)
integer           :: i

   do i = 1, NPar 
      !vals(i) = Npv_(i)*(MaxValue(i) - MinValue(i))  + MinValue(i)
      vals(i) = Npv_(i) * Ipv(i)
   enddo

end function Apv_
!-----------------------------------------------------------------------
subroutine EstimatePriors(C, InvC, error)
! To be called at the beginning of the run, 
! to estimate priors based on inital values for parameters
implicit none

real   , intent(OUT) :: C(NPar*(NPar+1)/2) 
real   , intent(OUT) :: InvC(NPar,NPar)
integer, intent(OUT) :: error
real                 :: work(NPar), U(NPar*(NPar+1)/2)
integer              :: m, q, nullty

! initialize the Prior Covariance matrix to zeros. 
  C(:) = 0d0
  do m = 1, NPar
! Assume a mean equal to the initial estimate for each parameter to be varied
!   nuPar and etaPar are also normalized values
     nuPar(m)     = Npv(m)

!   Assume a coefficient of variation
!     etaPar(m)    = CVprior * abs(nuPar(m))
     etaPar(m)    = 0.25*(MaxValue(m)-MinValue(m))/nuPar(m)
!   Set the variances on the diagonal of the Prior Covariance matrix
!   Here, assuming no cross-correlations
     C(m*(m+1)/2) = etaPar(m)**2  !Variance = sd**2
  enddo

! To output the parameter mean
  write(6,*) ' The prior parameter mean values: '
  write(6,1100) (nuPar(q), q = 1, Np2Vary)
! To output the parameter SD
  write(6,*) ' The prior SD of parameters: '
  write(6,1100) (etaPar(q), q = 1, Np2Vary)

! Invert the symmetric covariance matrix 
  call syminv(C, Np2Vary, U, work, nullty, error )

! Unpack the compressed inverse matrix U into the full (2-D) inverse covariance matrix
  call unpack(Np2Vary,U,InvC)

! To output the inverse prior covariance matrix
  write(6,*) ' The Inverse of the Prior Covariance Matrix: '
  do m = 1, Np2Vary
     write(6,1100) (InvC(m,q), q = 1, Np2Vary)
  end do

1100 format(<Np2Vary>(1pe9.1,1x))

! calculate the determinant of the prior covariance matrix (needed for the complete likelihood)

  call cholesky(C,Np2Vary,Np2Vary*(Np2Vary+1)/2,U,nullty,error)

! Just the products of the diagonals of the prior covariance matrix
  rtDetPcvm = 1d0
  do m = 1, Np2Vary
     rtDetPcvm = rtDetPcvm * U(m*(m+1)/2)
  end do

  write(6, *) ' The square root of the determinant of',      &
  ' the Prior Covariance Matrix = ', rtDetPcvm

END SUBROUTINE EstimatePriors
!-----------------------------------------------------------------------
subroutine SetSigmaLabels
implicit none
integer k

do k = 1, NDTYPE
   SigmaLabel(k) = trim('sigma')//DataLabel(k)
    SSqELabel(k) = trim('SSqE' )//DataLabel(k)
end do

end subroutine SetSigmaLabels
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
function CalcLogLike(ss, sig, param) result( lprob)
! Eq. 20 in Laine 2008 Ph.D thesis. This is the key subroutine!
implicit none
real, intent(in ) :: ss(NDTYPE)        ! sum of squared errors, for each data type
real, intent(in ) :: sig(NDTYPE)       ! sigmas for each data type
real, intent(in ) :: param(Np2Vary)    ! parameter values that can be varied
real              :: lprob             ! LogLikelihood

real              :: devpv(Np2Vary)
integer           :: k
  lprob = 0d0
  do k = 1, NDTYPE  !For different data types
! The line that is commented out includes the value of sigma for each data type. 
! With Gibbs sampling of the sigmas, the log probabilities with the newly proposed and current 
! parameter sets are compared, applying the same value of sigma for both. 
! Therefore, this routine should be called with the same current value of sigma, separately for
! the New and Current sum of square errors. 
!sig(k) : sigmas for each data type
!ss(k)  : sum of squared errors for each data type
!rt2pi  : square root of 2*pi

     lprob = lprob - NDPTS(k)*log(rt2Pi*sig(k)) - ss(k)/(2d0*sig(k)**2)
!!$       lprob = lprob - ss(k)/(2.0d0*sig(k)**2)
  enddo
!  write(6,*) 'lprob due to data-model mismatch: ', lprob
!!$ Including the effects of the priors
!!$ Assuming indepedent distributions: the kth parameter having mean = nuPar(k), sd = etaPar(k)
!!$    do k = 1, NPar
!!$       lprob = lprob - ((param(k) -nuPar(k))/(2.0d0*etaPar(k))**2) 
!!$    end do
!!$ Accounting for correlations between parameters in the prior
!!$ calculating the prior sum of squares using the inverse prior covariance matrix
    
    ! devpv: deviation of the parameter value from the mean
    devpv = param - nuPar  !param and nuPar must have the same length!!!

    ! rtDetPcvm: square root of the determinant of the prior covariance matrix
    lprob = lprob - log(rtDetPcvm*(rt2Pi**Np2Vary))             &
                  - sum( devpv*matmul(InvPriorCvm,devpv) )/2d0
!!     The later part equals to Eq. 23 in Laine thesis.
!!$    lprob = lprob - sum( devpv*matmul(InvPriorCvm,devpv) )/2.0d0

! for testing 
!!$    write(6, *) ' k:        param         devpv         nuPar       etaPar '
!!$    do k = 1, Np2Vary
!!$       write(6,1200) k, param(k), devpv(k), nuPar(k), etaPar(k)
!!$    END do
!!$    write(6, *) ' Prior contribution to lprob = ',sum( devpv*matmul(InvPriorCvm,devpv) )/2.0d0
!!$1200 format(1x,i4,2x,10(1pe12.3,2x) )

 return
end function CalcLogLike

!-----------------------------------------------------------------------
subroutine unpack(n,cA,fA)
implicit none
integer, intent(in)  :: n
real   , intent(in)  :: cA(n*(n+1)/2) ! compacted symmetric matrix A
real   , intent(out) :: fA(n,n)       ! fully expanded symmetric matrix A
integer              :: i, j, k

    k = 0
    do j = 1, n
      do i = 1, j - 1
        k = k + 1
        fA(i,j) = cA(k)
        fA(j,i) = cA(k)
      end do
      k = k + 1
      fA(j,j) = cA(k)
    end do

  return
end subroutine unpack
!-----------------------------------------------------------------------
END MODULE Interface_MOD
