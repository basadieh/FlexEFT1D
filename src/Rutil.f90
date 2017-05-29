module cholf90
  implicit none
contains
  subroutine cholesky ( a, n, nn, u, nullty, ifault )
!*****************************************************************************80
!
!! CHOLESKY computes the Cholesky factorization of a PDS matrix.
!
!  Discussion:
!
!    For a positive definite symmetric matrix A, the Cholesky factor U
!    is an upper triangular matrix such that A = U' * U.
!
!    This routine was originally named "CHOL", but that conflicted with
!    a built in MATLAB routine name.
!
!  Modified:
!
!    01 February 2008
!
!  Author:
!
!    Michael Healy
!    Modifications by AJ Miller.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Michael Healy,
!    Algorithm AS 6:
!    Triangular decomposition of a symmetric matrix,
!    Applied Statistics,
!    Volume 17, Number 2, 1968, pages 195-197.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A((N*(N+1))/2), a positive definite matrix 
!    stored by rows in lower triangular form as a one dimensional array, 
!    in the sequence
!    A(1,1),
!    A(2,1), A(2,2),
!    A(3,1), A(3,2), A(3,3), and so on.
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, integer NN, the dimension of the array used to store A, 
!    which should be at least (N*(N+1))/2.
!
!    Output, real ( kind = 8 ) U((N*(N+1))/2), an upper triangular matrix,
!    stored by columns, which is the Cholesky factor of A.  The program is
!    written in such a way that A and U can share storage.
!
!    Output, integer ( kind = 4 ) NULLTY, the rank deficiency of A.  If NULLTY is zero,
!    the matrix is judged to have full rank.
!
!    Output, integer ( kind = 4 ) IFAULT, an error indicator.
!    0, no error was detected;
!    1, if N < 1;
!    2, if A is not positive semi-definite.
!    3, NN < (N*(N+1))/2.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) ETA, should be set equal to the smallest positive
!    value such that 1.0 + ETA is calculated as being greater than 1.0 in the
!    accuracy being used.
!
  implicit none

  integer ( kind = 4 ), intent(in) :: n
  integer ( kind = 4 ), intent(in) :: nn

  real    ( kind = 8 ), intent(in) :: a(n)
  real    ( kind = 8 ), parameter :: eta = 1.0D-19
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ), intent(out) :: ifault
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ), intent(out) :: nullty
  real    ( kind = 8 ), intent(out) :: u(nn)
  real    ( kind = 8 ) w
  real    ( kind = 8 ) x

  ifault = 0
  nullty = 0

  if ( n <= 0 ) then
    ifault = 1
    return
  end if

  if ( nn < ( n * ( n + 1 ) ) / 2 ) then
    ifault = 3
    return
  end if

  j = 1
  k = 0
  ii = 0
!
!  Factorize column by column, ICOL = column number.
!
  do icol = 1, n

    ii = ii + icol
    x = eta * eta * a(ii)
    l = 0
    kk = 0
!
!  IROW = row number within column ICOL.
!
    do irow = 1, icol

      kk = kk + irow
      k = k + 1
      w = a(k)
      m = j

      do i = 1, irow - 1
        l = l + 1
        w = w - u(l) * u(m)
        m = m + 1
      end do

      l = l + 1

      if ( irow == icol ) then
        exit
      end if

      if ( u(l) /= 0.0D+00 ) then

        u(k) = w / u(l)

      else

        u(k) = 0.0D+00

        if ( abs ( x * a(k) ) < w * w ) then
          ifault = 2
          return
        end if

      end if

    end do
!
!  End of row, estimate relative accuracy of diagonal element.
!
    if ( abs ( w ) <= abs ( eta * a(k) ) ) then

      u(k) = 0.0D+00
      nullty = nullty + 1

    else

      if ( w < 0.0D+00 ) then
        ifault = 2
        return
      end if

      u(k) = sqrt ( w )

    end if

    j = j + icol

  end do

  return
end subroutine cholesky
subroutine subchl ( a, b, n, u, nullty, ifault, ndim, det )

!*****************************************************************************80
! 
!! SUBCHL computes the Cholesky factorization of a (subset of a) PDS matrix.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    FORTRAN77 version by Michael Healy, PR Freeman
!    FORTRAN90 version by  John Burkardt
!
!  Reference:
!
!    PR Freeman,
!    Remark AS R44:
!    A Remark on AS 6 and AS7: Triangular decomposition of a symmetric matrix
!    and Inversion of a positive semi-definite symmetric matrix,
!    Applied Statistics,
!    Volume 31, Number 3, 1982, pages 336-339.
!
!    Michael Healy,
!    Algorithm AS 6:
!    Triangular decomposition of a symmetric matrix,
!    Applied Statistics,
!    Volume 17, Number 2, 1968, pages 195-197.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A((M*(M+1))/2), a positive definite matrix 
!    stored by rows in lower triangular form as a one dimensional array, 
!    in the sequence
!    A(1,1),
!    A(2,1), A(2,2),
!    A(3,1), A(3,2), A(3,3), and so on.  
!    In the simplest case, M, the order of A, is equal to N.
!
!    Input, integer ( kind = 4 ) B(N), indicates the order in which the
!    rows and columns of A are to be used.  In the simplest case, 
!    B = (1,2,3...,N).
!
!    Input, integer ( kind = 4 ) N, the order of the matrix, that is, 
!    the matrix formed by using B to select N rows and columns of A.
!
!    Output, real ( kind = 8 ) U((N*(N+1))/2), an upper triangular matrix,
!    stored by columns, which is the Cholesky factor of A.  The program is
!    written in such a way that A and U can share storage.
!
!    Output, integer ( kind = 4 ) NULLTY, the rank deficiency of A.  
!    If NULLTY is zero, the matrix is judged to have full rank.
!
!    Output, integer ( kind = 4 ) IFAULT, an error indicator.
!    0, no error was detected;
!    1, if N < 1;
!    2, if A is not positive semi-definite.
!
!    Input, integer ( kind = 4 ) NDIM, the dimension of A and U, which might 
!    be presumed to be (N*(N+1))/2.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndim

  real    ( kind = 8 ) a(ndim)
  integer ( kind = 4 ) b(n)
  real    ( kind = 8 ) det
  real    ( kind = 8 ), parameter :: eta = 1.0D-09
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nullty
  real    ( kind = 8 ) u(ndim)
  real    ( kind = 8 ) w
  real    ( kind = 8 ) x

  ifault = 0
  nullty = 0
  det = 1.0D+00

  if ( n <= 0 ) then
    ifault = 1
    return
  end if

  ifault = 2
  j = 1
  k = 0

  do icol = 1, n

    ij = ( b(icol) * ( b(icol) - 1 ) ) / 2
    ii = ij + b(icol)
    x = eta * eta * a(ii)
    l = 0

    do irow = 1, icol

      kk = ( b(irow) * ( b(irow) + 1 ) ) / 2
      k = k + 1
      jj = ij + b(irow)
      w = a(jj)
      m = j

      do i = 1, irow - 1
        l = l + 1
        w = w - u(l) * u(m)
        m = m + 1
      end do

      l = l + 1

      if ( irow == icol ) then
        exit
      end if

      if ( u(l) /= 0.0D+00 ) then

        u(k) = w / u(l)

      else

        if ( abs ( x * a(kk) ) < w * w ) then
          ifault = 2
          return
        end if

        u(k) = 0.0D+00

      end if

    end do

    if ( abs ( eta * a(kk) ) <= abs ( w ) ) then

      if ( w < 0.0D+00 ) then
        ifault = 2
        return
      end if

      u(k) = sqrt ( w )

    else

      u(k) = 0.0D+00
      nullty = nullty + 1

    end if

    j = j + icol
    det = det * u(k) * u(k)

  end do

  return
end subroutine subchl
end module cholf90
!____________________________________________________________________________
! A C-program for MT19937: Real number version
!   genrand() generates one pseudorandom real number (double)
! which is uniformly distributed on [0,1]-interval, for each
! call. sgenrand(seed) set initial values to the working area
! of 624 words. Before genrand(), sgenrand(seed) must be
! called once. (seed is any 32-bit integer except for 0).
! Integer generator is obtained by modifying two lines.
!   Coded by Takuji Nishimura, considering the suggestions by
! Topher Cooper and Marc Rieffel in July-Aug. 1997.
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Library General Public
! License as published by the Free Software Foundation; either
! version 2 of the License, or (at your option) any later
! version.
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
! You should have received a copy of the GNU Library General
! Public License along with this library; if not, write to the
! Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
! 02111-1307  USA
!
! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
! When you use this, send an email to: matumoto@math.keio.ac.jp
! with an appropriate reference to your work.
!
!***********************************************************************
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
!
!   genrand()      -> double precision function grnd()
!   sgenrand(seed) -> subroutine sgrnd(seed)
!                     integer seed
!
! This program uses the following non-standard intrinsics.
!   ishft(i,n): If n>0, shifts bits in i by n positions to left.
!               If n<0, shifts bits in i by n positions to right.
!   iand (i,j): Performs logical AND on corresponding bits of i and j.
!   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
!
!***********************************************************************
! Fortran version rewritten as an F90 module and mt state saving and getting
! subroutines added by Richard Woloshyn. (rwww@triumf.ca). June 30, 1999

module mtmod
  implicit none
  ! Default seed
  integer, parameter :: defaultsd = 4357
  ! Period parameters
  integer, parameter :: N = 624, N1 = N + 1
  
  ! the array for the state vector
  integer, save, dimension(0:N-1) :: mt
  integer, save                   :: mti = N1
  
contains
  
  !Initialization subroutine
  subroutine sgrnd(seed)
    implicit none
!
!      setting initial seeds to mt[N] using
!      the generator Line 25 of Table 1 in
!      [KNUTH 1981, The Art of Computer Programming
!         Vol. 2 (2nd Ed.), pp102]
!
    integer, intent(in) :: seed

    mt(0) = iand(seed,-1)
    do mti=1,N-1
      mt(mti) = iand(69069 * mt(mti-1),-1)
    enddo
!
    return
  end subroutine sgrnd

!Random number generator
  real(8) function grnd()
    implicit integer(a-z)

! Period parameters
    integer, parameter :: M = 397, MATA  = -1727483681
!                                    constant vector a
    integer, parameter :: LMASK =  2147483647
!                                    least significant r bits
    integer, parameter :: UMASK = -LMASK - 1
!                                    most significant w-r bits
! Tempering parameters
    integer, parameter :: TMASKB= -1658038656, TMASKC= -272236544

    integer :: mag01(0:1)
    data mag01/0, MATA/
    save mag01
!                        mag01(x) = x * MATA for x=0,1

    TSHFTU(y)=ishft(y,-11)
    TSHFTS(y)=ishft(y,7)
    TSHFTT(y)=ishft(y,15)
    TSHFTL(y)=ishft(y,-18)

    if(mti.ge.N) then
!                       generate N words at one time
      if(mti.eq.N+1) then
!                            if sgrnd() has not been called,
        call sgrnd( defaultsd )
!                              a default initial seed is used
      endif

      do kk=0,N-M-1
          y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
          mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
      enddo
      do kk=N-M,N-2
          y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
          mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
      enddo
      y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
      mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
      mti = 0
    endif

    y=mt(mti)
    mti = mti + 1 
    y=ieor(y,TSHFTU(y))
    y=ieor(y,iand(TSHFTS(y),TMASKB))
    y=ieor(y,iand(TSHFTT(y),TMASKC))
    y=ieor(y,TSHFTL(y))

    if(y .lt. 0) then
      grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
    else
      grnd=dble(y)/(2.0d0**32-1.0d0)
    endif

    return
  end function grnd

 end module mtmod

!!$ program main
!!$! this main() outputs first 1000 generated numbers
!!$!
!!$    use mtmod
!!$    implicit none
!!$ 
!!$    integer, parameter      :: no=1000
!!$    real(8), dimension(0:7) :: r
!!$    integer j,k
!!$!    real(8) grnd
!!$!
!!$!      call sgrnd(4357)
!!$!                         any nonzero integer can be used as a seed
!!$    do j=0,no-1
!!$      r(mod(j,8))=grnd()
!!$      if(mod(j,8).eq.7) then
!!$        write(*,'(8(f9.6,'' ''))') (r(k),k=0,7)
!!$      else if(j.eq.no-1) then
!!$        write(*,'(8(f9.6,'' ''))') (r(k),k=0,mod(no-1,8))
!!$      endif
!!$    enddo
!!$
!!$ end program main
module gammaf90
  use mtmod
  implicit none

contains
 real(8) function gammadev(alpha,beta)  result(g)
!!$ Generates pseudo-random numbers from a Gamma distribution with parameters alpha and beta, 
!!$ such that the mean of the distribution, mu = alpha * beta. 
!!$ (Different notations are used for the second parameter, sometimes defined as lambda = 1/beta)
!!$ The algorithms in this routine are from "Random Number Generation and Monte Carlo Methods", by James E. Gentle
!!$ (2nd edition, Springer-Verlag, New York, 2003)
   implicit none
   real(8), intent(in) :: alpha, beta
   real(8) :: am1, b, u1, u2, v, t, x, y
   integer :: i, j

   am1 = alpha - 1.0d0
   g = -1

   if( alpha .gt. 1.0d0 ) then  
      do while( g .le. 0.0d0 ) 
         u1 = grnd()
         u2 = grnd()
         v = (alpha - 1.0d0/(6.0d0*alpha))*u1/(am1*u2)
         if( 2.0d0 * (u2 - 1.0d0) /am1 + v + 1.0d0/v .le. 2.0d0 ) then
            g = am1 * v
         else if( 2*log(u2)/am1 - log(v) + v .le. 1.0d0 ) then
            g = am1 * v
         end if
      end do
   else ! alpha <= 1
      t = 0.07 + 0.75*sqrt(1.0d0 - alpha)
      b = 1.0d0 + exp(-t)*alpha/t
      do while( g .le. 0.0d0 ) 
         u1 = grnd()
         u2 = grnd()
         v = b * u1
         if( v .le. 1 ) then
            x = t * v**(1.0d0/alpha)
            if( u2 .lt. (2.0d0-x)/(2.0d0+x) ) then
               g = x
            else if( u2 .lt. exp(-x) ) then
               g = x
            end if
         else
            x = - log(t*(b-v)/alpha)
            y = x/t
            if( u2*(alpha + y*(1.0d0-alpha)) .le. 1.0d0 ) then
               g = x
            else if( u2 .lt. y**am1 ) then
               g = x
            end if
         end if
      end do
   end if
   g = beta * g
   return
 end function gammadev
end module gammaf90
module gasdevf90
  implicit none

contains
  real(8) function gasdev() result(gval)
    use mtmod
    implicit none
    real(8) :: FAC, R, V1, V2, X
    real(8), save :: GSET
    integer, save :: ISET = 0

    IF (ISET.EQ.0) THEN
       R = 99 
       do while( R .ge. 1.0d0 )
          V1= 2.0d0*grnd() - 1.0d0
          V2= 2.0d0*grnd() - 1.0d0
          R = V1**2 + V2**2
       end do
       FAC = SQRT( -2.0d0*LOG(R)/R )
       GSET   = V1*FAC
       gval = V2*FAC
       ISET=1
    ELSE
       gval=GSET
       ISET=0
    ENDIF

    RETURN
  END function gasdev

end module gasdevf90
module mGf90
  use cholf90
  use gasdevf90
  implicit none

contains
 real(8) function multigauss(R,mu,n) result(G)
!!$ Generates a one dimensional array of length n containing multivariate Gaussian pseudo-random numbers
!!$ where R is the lower-triangular Cholesky factor of the covariance matrix of the desired distribution 
!!$ and mu is a one dimensional array of length n, containing the mean value for each random variable. 
!!$ R must be a one dimensional array, containing the lower-triangular matrix stored by row, starting from the 
!!$ uppermost, leftmost entry (first row, first column). 
   implicit none
   integer, intent(in) :: n
   real(8), intent(in) :: R(n*(n+1)/2)
   real(8), intent(in) :: mu(n)
   real(8) :: Nu(n)
   dimension :: G(n)
   integer :: i, j, nullty, cherror = 0
   
   if( n*(n+1)/2 .ne. size(R) ) then
      write(6,*) ' n*(n+1)/2 != size(R) in multiGauss '
      write(6,*) ' Cholesky factor matrix has size ', size(R) 
      write(6,*) ' but you are requesting a vector of size ',n
      stop
   else
      
!!$ generate an array of independent Gaussian variates
      do j = 1, n
         Nu(j) = gasdev()
      end do

!!$ start with the desired mean, and add the product of 
!!$ [the lower-triangular Cholesky factor of the covariance matrix] x [ the vector of iid Gaussian variates]
      do i = 1, n
         G(i) = mu(i)
         do j =  1,i
            G(i) = G(i) + R( i*(i-1)/2 + j ) * Nu(j)
         end do
      end do

   end if
   return
 end function multigauss
end module mGf90
module syminvf90
  use cholf90
contains
  subroutine syminv ( a, n, c, w, nullty, ifault )
    !*****************************************************************************80
    !
    !! SYMINV computes the inverse of a symmetric matrix.
    !
    !  Modified:
    !
    !    11 February 2008
    !
    !  Author:
    !
    !    FORTRAN77 version by Michael Healy
    !    FORTRAN90 version by John Burkardt
    !
    !  Reference:
    !
    !    Michael Healy,
    !    Algorithm AS 7:
    !    Inversion of a Positive Semi-Definite Symmetric Matrix,
    !    Applied Statistics,
    !    Volume 17, Number 2, 1968, pages 198-199.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A((N*(N+1))/2), a positive definite matrix stored 
    !    by rows in lower triangular form as a one dimensional array, in the sequence
    !    A(1,1),
    !    A(2,1), A(2,2),
    !    A(3,1), A(3,2), A(3,3), and so on.
    !
    !    Input, integer ( kind = 4 ) N, the order of A.
    !
    !    Output, real ( kind = 8 ) C((N*(N+1))/2), the inverse of A, or generalized
    !    inverse if A is singular, stored using the same storage scheme employed
    !    for A.  The program is written in such a way that A and U can share storage.
    !
    !    Workspace, real ( kind = 8 ) W(N).
    !
    !    Output, integer ( kind = 4 ) NULLTY, the rank deficiency of A.  If NULLTY is zero,
    !    the matrix is judged to have full rank.
    !
    !    Output, integer ( kind = 4 ) IFAULT, error indicator.
    !    0, no error detected.
    !    1, N < 1.
    !    2, A is not positive semi-definite.
    !
    implicit none
    
    integer ( kind = 4 ), intent(in) :: n
    
    real    ( kind = 8 ), intent(in)  :: a((n*(n+1))/2)
    real    ( kind = 8 ), intent(out) :: c((n*(n+1))/2)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) icol
    integer ( kind = 4 ), intent(out) :: ifault
    integer ( kind = 4 ) irow
    integer ( kind = 4 ) j
    integer ( kind = 4 ) jcol
    integer ( kind = 4 ) k
    integer ( kind = 4 ) l
    integer ( kind = 4 ) mdiag
    integer ( kind = 4 ) ndiag
    integer ( kind = 4 ) nn
    integer ( kind = 4 ) nrow
    integer ( kind = 4 ), intent(out) :: nullty
    real    ( kind = 8 ) w(n)
    real    ( kind = 8 ) x
    
    ifault = 0
    
    if ( n <= 0 ) then
       ifault = 1
       return
    end if
    
    nrow = n
    !
    !  Compute the Cholesky factorization of A.
    !  The result is stored in C.
    !
    nn = ( n * ( n + 1 ) ) / 2
    
    call cholesky ( a, n, nn, c, nullty, ifault )
    
    if ( ifault /= 0 ) then
       return
    end if
    !
    !  Invert C and form the product (Cinv)' * Cinv, where Cinv is the inverse
    !  of C, row by row starting with the last row.
    !  IROW = the row number, 
    !  NDIAG = location of last element in the row.
    !
    irow = nrow
    ndiag = nn
    
    do
       !
       !  Special case, zero diagonal element.
       !
       if ( c(ndiag) == 0.0D+00 ) then
          
          l = ndiag
          do j = irow, nrow
             c(l) = 0.0D+00
             l = l + j
          end do
          
       else
          
          l = ndiag
          do i = irow, nrow
             w(i) = c(l)
             l = l + i
          end do
          
          icol = nrow
          jcol = nn
          mdiag = nn
          
          do
             
             l = jcol
             
             if ( icol == irow ) then
                x = 1.0D+00 / w(irow)
             else
                x = 0.0D+00
             end if
             
             k = nrow
             
             do while ( irow < k )
                
                x = x - w(k) * c(l)
                k = k - 1
                l = l - 1
                
                if ( mdiag < l ) then
                   l = l - k + 1
                end if
                
             end do
             
             c(l) = x / w(irow)
             
             if ( icol <= irow ) then
                exit
             end if
             
             mdiag = mdiag - icol
             icol = icol - 1
             jcol = jcol - 1
             
          end do
          
       end if
       
       ndiag = ndiag - irow
       irow = irow - 1
       
       if ( irow <= 0 ) then
          exit
       end if
       
    end do
    
    return
  end subroutine syminv
end module syminvf90
