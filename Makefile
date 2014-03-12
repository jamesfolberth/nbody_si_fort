
.PHONY: all build clean force

FORT=gfortran

FORTOPTS = -Wall -fbounds-check -funroll-loops -finline -ffast-math -march=native -mtune=native -O2

LDLIBS= -lsatlas

build:
	make -C src build


# phony target to force make to build targets that depend on force
force:

clean:
	rm -f main
	rm -rf *.o *.mod *.out *.a

