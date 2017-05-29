#!/bin/bash



export RUNDIR=~/Working/FlexEFT1D2/S1
cd $RUNDIR

/bin/cp ../bio_rhs.f90 $RUNDIR
/bin/cp ../FlexEFT.f90 $RUNDIR

ifort -o flex bio_rhs.f90 FlexEFT.f90
