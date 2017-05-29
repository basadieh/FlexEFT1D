subroutine MCMC_adapt(Nmcmc, intv, ID_)
USE IO_MOD
! Major MCMC subroutine
implicit none
integer, intent(in)       ::  Nmcmc   ! The time steps for random MCMC chains

logical                   ::  BURNIN = .FALSE.
! The number of steps for adapting the PCVM
! if interval >= N (i.e. EnsLen <= 1), does not adapt at all
integer, intent(in)       ::  intv, ID_ 


real   :: SS1(NDTYPE), SS2(NDTYPE), SS3(NDTYPE), ALPHA12, alpha13
real   :: NewLogLike2
real   :: PCVM2(NPar*(NPar+1)/2)
real   :: subppro2(NPar)
real   :: AMaccR, DRaccR
! Output also includes:
! 1) Params, CVM, PCVM, SSqE, sigma, loglikehood at each step
! 2) Best of the above quantities so far
integer :: i,j, k, EnsLen

! Downscaling factor for Delayed Rejection MCMC
real,    parameter        :: DRScale = 1d-2


! The number of times of adaptation (i.e. number of ensembles)
EnsLen = INT(Nmcmc/intv)

If (EnsLen .LE. 1) Then
   BURNIN   = .TRUE.
   EnsLen   = 1
Endif
! Acceptance counts for Adaptive MCMC
AMacc = 0

! Acceptance counts for Delayed Rejection MCMC
DRacc = 0

BestLogLike = -1d12  ! A very large, negative number for very low probability
MeanLogLike = 0d0

DO i = 1, EnsLen

!!$  write output to Ensemble files for simulated values
   call model_ensout(eofint,jrun,ID_,subpcurr, dumE)

!  write best results into the best file
   call write_best(ID_)

   subcycle: do j = 1, intv
     
               jrun = jrun + 1
           subpmean =((real(jrun)-1d0)*subpmean + subpcurr)    /real(jrun)
           sigmamean=((real(jrun)-1d0)*sigmamean+    sigma)    /real(jrun)
         MeanLogLike=((real(jrun)-1d0)*MeanLogLike+CurrLogLike)/real(jrun)

      if (verbose) then
         if (BURNIN) then
            write(6,*) 'Burn-in ', jrun, 'th iteration:'
         else
            write(6,*) 'Main loop ', jrun, 'th iteration:'
         endif
      endif

      ! write current parameters and sigma into output file
      Call write_ensemble

! Calculate the probability of accepting the new PAR

  ! Calcuate the old SSqE
    SS1   = CurrSSqE

  ! Calcuate the new SSqE
    SS2   = CalSSQE(subppro)

  ! subroutine MCMC_alpha(newSS, newpar, oldSS, oldpar, newL, oldL, alpha) 
    call MCMC_alpha(SS2, subppro, SS1, subpcurr,  &
                    NewLogLike, CurrLogLike, ALPHA12)

! The first step to evaluate whether to accept new parameters
    IF( GRND() .LT. ALPHA12 ) THEN
       !accept the first move (Y1)
       !acc is a tally of the acceptance rate, the rest are parameters
      ! Counting of accepted times
      AMacc       = AMacc+1
      CurrLogLike = NewLogLike
      CurrSSqE    = SS2
      subpcurr    = subppro

   
    ELSE  ! REJECT THE FIRST MOVE (Y1)
    ! DR:
    ! Propose a second move (Y2) based on scaled proposal covariance matrix
    ! and current position
      if (DRAM) then

          PCVM2      = DRScale*PCVM
          subppro2   = NewPAR(PCVM2, subpcurr)

    !      Calcuate the new SSqE
           SS3       = CalSSQE(subppro2)

    ! Calculate the acceptance probability of the second move
!subroutine MCMC_DR_alpha13(C, &
!                           par1, ss1, logLike1, &
!                           par2, ss2, logLike2, &
!                           par3, ss3, logLike3, alpha13)

           call MCMC_DR_alpha13(PCVM,                        &
                           subpcurr, SS1, CurrLogLike,  &
                           subppro,  SS2,  NewLogLike,  &
                           subppro2, SS3,  NewLogLike2, alpha13)

          ! Judge whether the second move should be accepted or not
            if( GRND() .lt. alpha13 ) then
                !accept the second move (Y1)

               DRacc       = DRacc+1
               CurrLogLike = NewLogLike2
               CurrSSqE    = SS3
               subpcurr    = subppro2

            endif

      endif

    ENDIF  !<= THE ENDIF OF ACCEPTING NEW PARAMTERS.

    !keep track of the 'best' run
    if( CurrLogLike .gt. (BestLogLike+1d-10) ) then
       write(6,*) 'Get best Log-Likelihood!'
       BestLogLike = NewLogLike     
       BestSSqE    = CurrSSqE
       sigmabest   = sigma
       subpbest    = subpcurr
    endif


    ! Update the Covariance matrix and the Parameter mean
    ! subroutine UpdateCVM(OldCVM, Oldmean, n, par, newmean, NewCVM)
    call UpdateCVM(Cvm, subpcurrmean, DBLE(jrun), subpcurr,      &
                   subpcurrmean, Cvm)


    ! Propose newpar based on old par
    subppro = NewPAR(Pcvm, subpcurr)

    ! Update sigma:
    sigma = newsigma(CurrSSqE)

   enddo subcycle

! Calculate new PCVM based on CVM

!!$  Set the Proposal covariance matrix, by scaling the Parameter Covariance matrix
     Pcvm = Cvm*Spcvm/NPar
!!$ Add a small term to the diagonal entries, so that the matrix will not be singular. 
     do k = 1, NPar
       Pcvm(k*(k+1)/2)=Pcvm(k*(k+1)/2) + CvEpsilon*Spcvm/NPar
     enddo

ENDDO

AMaccR = real(AMacc)/real(Nmcmc)
DRaccR = real(DRacc)/real(Nmcmc)

If (BURNIN) Then
   write(6,101) 'The acceptance rate of 1st move during Burn-in is: ', AMaccR
   write(6,101) 'The acceptance rate of 2nd move during Burn-in is: ', DRaccR
Else
   write(6,101) 'The acceptance rate of 1st move during the major loop is: ', AMaccR
   write(6,101) 'The acceptance rate of 2nd move during the major loop is: ', DRaccR
Endif
101 Format(A100, 1F12.3)
End subroutine MCMC_adapt

