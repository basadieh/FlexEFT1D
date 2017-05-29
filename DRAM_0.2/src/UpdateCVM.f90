subroutine UpdateCVM(OldCVM, Oldmean, n, par, newmean, NewCVM)
Use Interface_MOD, only: NPar, PriorCVM
implicit none
real, intent(in)   :: OldCVM(NPar*(NPar+1)/2)
real, intent(in)   :: Oldmean(NPar)
real, intent(in)   :: n  !current step as a weight
real, intent(in)   :: par(NPar)
real, intent(out)  :: newmean(NPar)
real, intent(out)  :: NewCVM(NPar*(NPar+1)/2)

integer            :: row,col

IF (n .GE. 2d0) THEN
   ! Update the mean params
   
   newmean = ((n-1d0)*Oldmean + par)/n
   
   ! Update the correlation covariance matrix (but not Pcvm):    
   do row = 1, NPar
      do col = 1, row 
        NewCVM(row*(row-1)/2 +col) =                  & 
            ((n-1d0)/n)*OldCvm(row*(row-1)/2 +col)    &
          + (1d0/(n-1d0))                             &
          * (par(row)-newmean(row))                   &
          * (par(col)-newmean(col))
      enddo
   enddo
ELSE
   newmean   = par
   NewCVM(:) = PriorCvm(:)
ENDIF
End subroutine

