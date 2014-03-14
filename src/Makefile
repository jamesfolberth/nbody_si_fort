
.PHONY: all build clean force
SHELL=/usr/bin/bash

FORT=gfortran

FORTOPTS = -Wall -fbounds-check -funroll-loops -finline -ffast-math -march=native -mtune=native -O2

PWD=$(shell pwd)
HDF5_HOME=$(PWD)/../hdf5
INCHDF= -I$(HDF5_HOME)/include
LIBHDF= -L$(HDF5_HOME)/lib -lhdf5_fortran

LDLIBS= -lsatlas
LDADDPATH= -Wl,-rpath=$(HDF5_HOME)/lib # add this dir to the runtime search path for ld


all: orbit h5_crtdat

utils.o: utils.f90
	$(FORT) -c utils.f90 $(FORTOPTS) $(INCHDF)

orbital_data.o: orbital_data.f90
	$(FORT) -c orbital_data.f90 $(FORTOPTS) $(INCHDF)

coords.o: coords.f90
	$(FORT) -c coords.f90 $(FORTOPTS) $(INCHDF)

kepler.o: kepler.f90
	$(FORT) -c kepler.f90 $(FORTOPTS) $(INCHDF)

orbit: orbit.f90 utils.o orbital_data.o coords.o kepler.o
	$(FORT) -c orbit.f90 $(FORTOPTS) $(INCHDF)
	$(FORT) -o orbit $(FORTOPTS) $(LIBHDF) $(LDLIBS) $(LDADDPATH) utils.o orbital_data.o coords.o kepler.o orbit.o 
	./orbit

# use ldd h5_crtdat to check that the proper lib is going to load
h5_crtdat: h5_crtdat.f90
	$(FORT) -c h5_crtdat.f90 $(FORTOPTS) $(INCHDF) $(LIBHDF)
	$(FORT) -o h5_crtdat $(FORTOPTS) $(LIBHDF) $(LDLIBS) $(LDADDPATH) h5_crtdat.o

run_h5_crtdat: h5_crtdat
	./h5_crtdat

# phony target to force make to build targets that depend on force
force:

clean:
	rm -rf *.o *.mod *.out *.a

