!-----------------------------------------------------------------------
subroutine transform(N,x,xmin,xmax,convert,y)
use Interface_MOD, only : NO
implicit none

integer, intent(in)       :: N

! Input vector can be original values or need to be converted back
real(8), intent(in)       :: x(N)

! The minimal allowed positive values in x, to avoid zero or negative values
! These values are normal values before transformation
real(8), intent(in)       :: xmin,xmax

! Whether to lognorm or convert back to normal values
integer, intent(in)       :: convert

real(8), intent(out)      :: y(N)

real(8) :: max_y,min_y

integer :: i

! Normalization occurs after sqrt(sqrt(x)))  !
min_y = sqrt(sqrt(xmin))
max_y = sqrt(sqrt(xmax))

if (convert .eq. NO) then
   do i = 1, N
      y(i) = sqrt(sqrt(x(i)))
      y(i) = (y(i)-min_y)/(max_y-min_y)
   enddo 
else
   do i = 1, N
      y(i) = x(i)*(max_y-min_y) + min_y
      y(i) = y(i)**4
   enddo 
endif
end subroutine transform

