CFT=ifort
.SUFFIXES: .o .i .f90 .f
.f90.o:
	$(CFT) -c -O3 -r8 -i4 $*.f90 -o $*.o
.f.o:
	$(CFT) -c -O3 -r8 -i4 $*.f -o $*.o

MAIN = DRAM
#Source file
SRCS = dsymv.f lsame.f xerbla.f Readcsv.f90\
	cholesky.f90 gridinterp.f90 time_interp.f90\
	mtfort90.f90 gammadev.f90 gasdev.f90 multiGauss.f90\
	bio_MOD.f90 MOD_1D.f90 \
	matmuls.f90\
	syminv.f90 Interface_MOD.f90 transform.f90 UpdateCVM.f90\
	SUB_MOD.f90 Main.f90

OBJS =$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SRCS)))
#Executable file
$(MAIN): $(OBJS)
	$(CFT) -O3 -r8 -i4 -o $(MAIN) $(OBJS)

