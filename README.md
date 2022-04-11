# simdsort

A fast quicksort function with AVX2,  AVX512 and SVE  support

you can try this with

```
  make testsort  #(C version)
```
or
```
  make cpptestsort  #(C++ version)
```
or
```
  make kvtestsort  #(C++ key-value pair version)
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

On Fugaku, the same test resulted in:
```
qsort:  0.00151181
C++ std::sort:   0.000728001
blocksort:  0.00029474
sort passed for n=8002
```

## Usage as header-only library

### Pure C:

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

Current Makefile check the architecture using /usr/bin/arch. So try to
use native compiler on Fugaku (not the cross compiler). 

## Note

This source is in pure C and no care for name collision is done.


### C++

```
   #include "PATH/OF/THIS/FILE/simdsort.hpp"
```

Provide:
```
void simd_sort(int64_t * r, int n);
void simd_sort(uint64_t* hi,uint64_t* lo, uint64_t* index,  int n);
```
in namespace SIMDSortLib

The version for single int64_t array is identical as Pure C version.
r is the pointer to the array of int64_t to be sorted
n is the number of elements

The version which takes three arrays sorts these three arrays
regarding three words hi[i], lo[i], index[i] forming single 192-bit
unsigned integer. 

## Algorithm

Essentially the same as in
Fast Quicksort Implementation Using AVX Instructions,
Gueron and  Krasnov 2015, 
10.1093/comjnl/bxv063

In the AVX2, AVX512 and SVE  versions, bitonic sort is called for
n<=8, n<=16 and n<=16, respectively.  (192-bit version works only with AVX512)

## Changelog

### Mar 27, 2022

* AVX2 version uses  bitonic sort for n<=8

### Mar 29, 2022

* AVX512 version uses  bitonic sort for n<=16

### Mar 30, 2022

* SVE version uses  bitonic sort for n<=16

### Apr 8, 2022

* 192-bit uint version and C++ interface added

### Apr 10, 2022

* 192-bit uint version use AVX512/SVE and bitonic sort
* Makefile uses /usr/bin/arch to determine which compiler and flags to be used
