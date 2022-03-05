#% : %.cpp
#	g++ $(CPPFLAGS) -o $@ $<
% : %.c
	gcc $(CFLAGS) -o $@ $<
CFLAGS = -g  -O3 -mavx2 -ftree-vectorize   -fopt-info-vec
CPPFLAGS = -ffast-math    -ftree-vectorize  -fopt-info-vec-optimized=vector.txt



testsort: simdsorttest
	echo Test ssamplesort for n=2 to 8000. This will take a while...
	./simdsorttest 2 8000
	./simdsorttest 8001 8003 1
simdsorttest: simdsorttest.c simdsort.h
