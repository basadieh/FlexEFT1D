        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 29 14:19:28 2016
        MODULE GRIDINTERPOL__genmod
          INTERFACE 
            SUBROUTINE GRIDINTERPOL(N,COLS,OBS_Z,OBS_PROF,NLEV_,MODEL_Z,&
     &MODEL_PROF)
              INTEGER(KIND=4), INTENT(IN) :: NLEV_
              INTEGER(KIND=4), INTENT(IN) :: COLS
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: OBS_Z(N)
              REAL(KIND=8), INTENT(IN) :: OBS_PROF(N,COLS)
              REAL(KIND=8), INTENT(IN) :: MODEL_Z(NLEV_)
              REAL(KIND=8), INTENT(OUT) :: MODEL_PROF(NLEV_,COLS)
            END SUBROUTINE GRIDINTERPOL
          END INTERFACE 
        END MODULE GRIDINTERPOL__genmod
