
.PHONY: all build clean force
SHELL=/usr/bin/bash

all: orbit lyapunov

orbit: force
	make -C src orbit

lyapunov: force
	make -C src lyapunov

# Build odepack.so (F77)
odepack:
	make -C odepack

test_dset:
	make -C src h5_crtdat
	make -C src run_h5_crtdat

# phony target to force make to build targets that depend on force
force:

clean:
	make -C src clean
	

