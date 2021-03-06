
.PHONY: all build clean force
SHELL=/usr/bin/bash

FORT=gfortran

#FORTOPTS = -Wall -fbounds-check -funroll-loops -finline -ffast-math -march=native -mtune=native -O2
FORTOPTS = -Wall -Wextra -march=native -mtune=native -funroll-loops -finline -finline-functions -O2
#FORTOPTS = -Wall -Wextra -march=native -mtune=native -funroll-loops -finline -finline-functions -O0 -g -Q

PWD=$(shell pwd)
HDF5_HOME=$(PWD)/../hdf5
INCHDF= -I$(HDF5_HOME)/include
LIBHDF= -L$(HDF5_HOME)/lib -lhdf5_fortran

LDLIBS= -lsatlas
LDADDPATH= -Wl,-rpath=$(HDF5_HOME)/lib # add this dir to the runtime search path for ld

# odepack library
ODE_HOME=$(PWD)/../odepack
ODELIBS= -L$(ODE_HOME) -lodepack
ODEADDPATH= -Wl,-rpath=$(ODE_HOME)

all: orbit lyapunov orbit_lin_test

# Modules
coords.o: coords.f90 utils.o orbital_data.o
	$(FORT) -c coords.f90 $(FORTOPTS)

kepler.o: kepler.f90 utils.o
	$(FORT) -c kepler.f90 $(FORTOPTS)

orbital_data.o: orbital_data.f90 utils.o
	$(FORT) -c orbital_data.f90 $(FORTOPTS)

symplectic.o: symplectic.f90 utils.o
	$(FORT) -c symplectic.f90 $(FORTOPTS)

utils.o: utils.f90
	$(FORT) -c utils.f90 $(FORTOPTS) $(INCHDF)

# Main programs
orbit: orbit.f90 coords.o kepler.o orbital_data.o symplectic.o utils.o
	$(FORT) -c orbit.f90 $(FORTOPTS) $(INCHDF)
	$(FORT) -flto -o orbit $(FORTOPTS) $(LIBHDF) $(LDLIBS) $(LDADDPATH) utils.o orbital_data.o coords.o kepler.o symplectic.o orbit.o 
	nice -n0 ./orbit

lyapunov: lyapunov.f90 coords.o kepler.o orbital_data.o symplectic.o utils.o force
	$(FORT) -c lyapunov.f90 $(FORTOPTS) $(INCHDF)
	$(FORT) -flto -o lyapunov $(FORTOPTS) $(LIBHDF) $(LDLIBS) \
		$(LDADDPATH) $(ODELIBS) $(ODEADDPATH) \
		utils.o orbital_data.o coords.o kepler.o symplectic.o lyapunov.o 
	nice -n0 ./lyapunov

orbit_lin_test: orbit_lin_test.f90 coords.o kepler.o orbital_data.o symplectic.o utils.o
	$(FORT) -c orbit_lin_test.f90 $(FORTOPTS) $(INCHDF)
	$(FORT) -flto -o orbit_lin_test $(FORTOPTS) $(LIBHDF) $(LDLIBS) \
		$(LDADDPATH) $(ODELIBS) $(ODEADDPATH) \
		utils.o orbital_data.o coords.o kepler.o symplectic.o orbit_lin_test.o 
	#nice -n0 ./orbit_lin_test



## use ldd h5_crtdat to check that the proper lib is going to load
#h5_crtdat: h5_crtdat.f90
#	$(FORT) -c h5_crtdat.f90 $(FORTOPTS) $(INCHDF) $(LIBHDF)
#	$(FORT) -o h5_crtdat $(FORTOPTS) $(LIBHDF) $(LDLIBS) $(LDADDPATH) h5_crtdat.o
#
#run_h5_crtdat: h5_crtdat
#	./h5_crtdat

# phony target to force make to build targets that depend on force
force:

clean:
	rm -f h5_crtdat orbit lyapunov
	rm -f *.o *.mod *.out *.a

