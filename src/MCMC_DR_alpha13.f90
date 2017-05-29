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

use sub_mod, only: Invmat, NPar, CalcLogLike, sigma, NDTYPE
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
real             :: cff1(NPar), cff2(NPar)
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

  call matmuls(iC, (par3-par2), cff1) 
  call matmuls(iC, (par1-par2), cff2) 

!  q1 = -0.5d0*( &
!       sum(  matmuls(iC,(par3-par2)) * (par3-par2) ) - &
!       sum(  matmuls(iC,(par1-par2)) * (par1-par2) ) )

 q1 = -0.5d0*( &
      sum(  cff1 * (par3-par2) ) - &
      sum(  cff2 * (par1-par2) ) )


  ! alpha2(x,y1,y2)
  alpha13 = min(1d0, exp(l2+q1)*(1d0 - alpha32)/(1d0 - alpha12))

end subroutine MCMC_DR_alpha13
!-----------------------------------------------------------------------

