#% : %.cpp
#	g++ $(CPPFLAGS) -o $@ $<
% : %.c
	gcc $(CFLAGS) -o $@ $<
CFLAGS = -g -DAVX512 -O3 -mavx512f -ftree-vectorize   -fopt-info-vec # -pg
#CFLAGS = -g -DAVX2 -O3 -mavx2 -ftree-vectorize   -fopt-info-vec # -pg
CC  = gcc
#CFLAGS = -g -Nclang  -Ofast -msve-vector-bits=512 -mcpu=a64fx+sve
#CC = FCCpx

simdsorttest: simdsorttest.c simdsort.h Makefile sort.cpp bitonic8.h bitonic16.h
	$(CC) $(CFLAGS) -o simdsorttest simdsorttest.c  sort.cpp
testsort: simdsorttest
	echo Test ssamplesort for n=2 to 8000. This will take a while...
	./simdsorttest 2 8000
	./simdsorttest 8001 8003 1
bitonic16.h: bitonic.rb
	ruby bitonic.rb AVX512 3 > bitonic16.h
