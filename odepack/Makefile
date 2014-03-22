
.PHONY: all clean
SHELL=/usr/bin/bash

FORT=gfortran

FORTOPTS = -Wall -Wextra -march=native -mtune=native -funroll-loops -finline -finline-functions -Ofast
#FORTOPTS = -g -O0


all: libodepack.so

libodepack.so: opkdmain.f opkda1.f opkda2.f
	$(FORT) $(FORTOPTS) -shared opkdmain.f opkda1.f opkda2.f \
		-o libodepack.so -fPIC

clean:
	rm -f libodepack.so
	
