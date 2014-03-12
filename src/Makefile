
.PHONY: all build clean force

FORT=gfortran

FORTOPTS = -Wall -fbounds-check -funroll-loops -finline -ffast-math -march=native -mtune=native -O2

HDF5_HOME=../hdf5
INCHDF= -I$(HDF5_HOME)/include
LIBHDF= -L$(HDF5_HOME)/lib -lhdf5_fortran

LDLIBS= -lsatlas
LDADDPATH= -Wl,-rpath=$(HDF5_HOME)/lib # add this dir to the runtime search path for ld


all: build

build: main.f90
	$(FORT) -c main.f90 $(FORTOPTS) $(INCHDF) $(LIBHDF)
	$(FORT) -o main $(FORTOPTS) $(LIBHDF) $(LDLIBS) $(LDADDPATH) main.o

h5_crtdat: h5_crtdat.f90 force
	$(FORT) -c h5_crtdat.f90 $(FORTOPTS) $(INCHDF) $(LIBHDF)
	$(FORT) -o h5_crtdat $(FORTOPTS) $(LIBHDF) $(LDLIBS) $(LDADDPATH) h5_crtdat.o # use ldd h5_crtdat to check that the proper lib is going to load


# phony target to force make to build targets that depend on force
force:

clean:
	rm -f main h5_crtdat
	rm -rf *.o *.mod *.out *.a

