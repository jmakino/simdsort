% : %.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<
% : %.c
	$(CC) $(CFLAGS) -o $@ $<

ARCH = $(shell /usr/bin/arch)
ifeq ($(ARCH),x86_64)
  CFLAGS = -g -DAVX512 -O3 -mavx512f -ftree-vectorize   -fopt-info-vec # -pg
  CC  = gcc
  CPP = g++
endif
ifeq ($(ARCH),aarch64)
  CFLAGS = -g -DSVE -Nclang  -Ofast -mcpu=a64fx+sve # -msve-vector-bits=512 
  CC  = fcc
  CPP = FCC
endif
CPPFLAGS = $(CFLAGS) 


EXES = simdsorttest sorttest kvsorttest
simdsorttest: simdsorttest.c simdsort.h Makefile sort.cpp bitonic8.h bitonic16.h  bitonic16sve.h
	$(CC) $(CFLAGS) -o simdsorttest simdsorttest.c  sort.cpp
testsort: simdsorttest
	echo Test simdsort for n=2 to 8000. This will take a while...
	./simdsorttest 2 8000
	./simdsorttest 8000 8003 1
cpptestsort: sorttest
	echo Test simdsort for n=2 to 8000. This will take a while...
	./sorttest 2 8000
	./sorttest 8000 8003 1
kvtestsort: kvsorttest
	echo Test kvsimdsort for n=2 to 8000. This will take a while...
	./kvsorttest 2 8000
	./kvsorttest 10000 10003 1
sorttest: sorttest.cpp simdsort.hpp
kvsorttest: kvsorttest.cpp simdsort.hpp bitonic16kv.h  bitonic16.h bitonic16svekv.h  bitonic16sve.h
bitonic16.h: bitonic.rb
	ruby bitonic.rb AVX512 3 > bitonic16.h
bitonic16kv.h: bitonic.rb
	ruby bitonic.rb AVX512KV 3 > bitonic16kv.h
bitonic16sve.h: bitonic.rb
	ruby bitonic.rb SVE 3 > bitonic16sve.h
bitonic16svekv.h: bitonic.rb
	ruby bitonic.rb SVEKV 3 > bitonic16svekv.h

bitonic4.h: bitonic.rb
	ruby bitonic.rb 4 3 > bitonic4.h
clean:
	rm $(EXES)

