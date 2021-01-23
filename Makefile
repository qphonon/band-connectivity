all : band-connectivity

#F90=ifort -fPIC -ftrapuv -warn all -CB
F90=gfortran -fPIC -ffree-line-length-0

#lapack=-L /install-linux-ifort/lapack-3.8.0 -llapack -lrefblas
#lapack=/usr/lib/x86_64-linux-gnu/lapack/liblapack.a /usr/lib/x86_64-linux-gnu/blas/libblas.a 
lapack=-L /install-linux-gfortran/lapack-3.8.0 -llapack -lrefblas


MYPROG=.

types_and_constants.o: $(MYPROG)/types_and_constants.f90
	$(F90) -c types_and_constants.f90

derived_constants.o: types_and_constants.o $(MYPROG)/derived_constants.f90
	$(F90) -c derived_constants.f90
	
sgsym.o: types_and_constants.o $(MYPROG)/sgsym.f90
	$(F90) -c sgsym.f90

commod.o: types_and_constants.o derived_constants.o sgsym.o $(MYPROG)/commod.f90
	$(F90) -c commod.f90

qmod.o : types_and_constants.o derived_constants.o sgsym.o commod.o 
	$(F90) -c qmod.f90
	
band-connectivity : qmod.o
	$(F90) -o $@ types_and_constants.o sgsym.o commod.o qmod.o band-connectivity.f90 $(lapack) 

clean:
	- rm fm-forces *.mod types_and_constants.o derived_constants.o sgsym.o commod.o
	- rm -rf tbEIGENVAL arm-chair simple eigenmode-*.axsf ampl-eigenmode-*.dat  vc-eigenmode-*-movie.axsf movie.axsf unformatted-phonon-dos.dat

