# simdsort

A fast quicksort function with AVX2,  AVX512 and SVE  support

you can try this with

```
  make testsort
```
The last part of the  output would look like
```
qsort:  0.000468138
C++ std::sort:   0.000298972
blocksort:  9.8489e-05
sort passed for n=8002
```
The numbers are elapsed time to sort 8002 64-bit integer numbers,
with system qsort(3), C++ std::sort  and simd quicksort.
This result is ontained  on Intel® Core™ i7-1065G7 with gcc version
7.5.0 (Ubuntu 7.5.0-3ubuntu1~18.04).   AVX2 version is used.

On Fugaku, the same test (slighyly older version) resulted in:
```
qsort:  0.00147074
sort:   0.000958731
blocksort:  0.000527481
sort passed for n=8002
```

## Usage as header-only library

```
   #include "PATH/OF/THIS/FILE/simdsort.h"
```



Provide:
```
void simd_sort_int64( int64_t * r, int n);
```


r is the pointer to the array of int64_t to be sorted
n is the number of elements

Use  -DAVX2, -DAVX512 and -DSVE   to utilize AVX2, AVX512 and SVE

For SVE, select CC = FCCpx in Makefile

## Note

This source is in pure C and no care for name collision is done.

## Algorithm

Essentially the same as in
Fast Quicksort Implementation Using AVX Instructions,
Gueron and  Krasnov 2015, 
10.1093/comjnl/bxv063

In the AVX2 version, bitonic sort is called for n<= 8.

## Changelog

### Mar 27, 2022

* AVX2 version uses  bitonic sort for n<=8

### Mar 29, 2022

* AVX512 version uses  bitonic sort for n<=16

